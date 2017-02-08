/* #ifndef DASH_ENABLE_LOGGING
#define DASH_ENABLE_LOGGING
#endif */

#include <libdash.h>
#include <chrono>
#include <thread>
#include "../Terminal_Color.h"

#define DEBUG
#define SLEEP_TIME            30 //that's the sleep time before DEBUG IO
#define MAX_KEY               99
#define MIN_NUM_ELEM_PER_UNIT 10


using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::pair;
using std::max;
using std::make_pair;
using std::this_thread::sleep_for;

using uint  = unsigned  int;
using POI_T = int;  //this type musst be signed!

template<typename T = POI_T>
inline void outer(         
                           uint   nelts ,
            vector< pair<T,T> > & points,
      dash::NArray< double, 2 > & matOut,
            dash::Array<double> & vec,
                            int   myid
    );

template< typename T >
void print2d( T& mat ) {
  for( int i = 0; i < mat.extent(0); i++ ) {
    for( int j = 0; j < mat.extent(1); j++ ) {
      cout << std::setw(3) << static_cast<uint>( mat(i,j) )<< " ";
    }
    cout << endl;
  }
}

template<typename T = POI_T>
inline void read_vector_of_points( uint nelts, vector< pair<T,T> > & points ) {
  for( uint i = 0; i < nelts; i++ ) {
    cin >> points[i].first >> points[i].second;
  }
}

template<typename T = POI_T>
inline void broadcastInputToUnits( uint nelts, vector< pair<T,T> > & points ) {
  dash::team_unit_t TeamUnit0ID = dash::Team::All().myid( );
  TeamUnit0ID.id = 0;
  dart_ret_t ret = dart_bcast(
                      static_cast<void*>( points.data( ) ),  // buf 
                      points.size( ) * sizeof(pair<T,T>)  ,  // nelts
                      DART_TYPE_BYTE                      ,  // dtype
                      TeamUnit0ID                         ,  // root
                      dash::Team::All().dart_id( )           // team
                   );
  if( DART_OK != ret ) cout << "An error while BCAST has occured!" << endl; 
}

int main( int argc, char* argv[] )
{  
  dash::init( &argc,&argv );
  int myid = static_cast<int>( dash::Team::GlobalUnitID( ).id );
  
  
  if( argc != 2 ){
    if( 0 == myid ){ cout << "1 Parameter expected!"           << endl
                          << "Usage: cowichan_outer nElements" << endl
                          << "Then enter the points."          << endl;
    }
    dash::finalize( );
    return 0;
  }
  
  uint nelts = static_cast<uint>( atoi( argv[1] ) );
  
  vector< pair<POI_T, POI_T> > points( nelts );
  dash::NArray< double, 2 >    matOut( nelts, nelts );
  dash::Array < double >       vec( nelts );
  
  //read input points on unit 0
  if( 0 == myid ) read_vector_of_points( nelts, points );
  
  //broadcast input from unit 0
  broadcastInputToUnits( nelts, points );

  outer( nelts, points, matOut, vec, myid );

  dash::finalize( );
}

inline double sqr(double const x) {
  return x * x;
}

template<typename T = POI_T>
inline double distance(const pair<T, T>& x, const pair<T, T>& y) {
  return sqrt(sqr(x.first - y.first) + sqr(x.second - y.second));
}

template<typename T = POI_T>
inline void outer(         
                           uint   nelts ,
            vector< pair<T,T> > & points,
      dash::NArray< double, 2 > & matOut,
            dash::Array<double> & vec,
                            int   myid
    ){

  dash::Team & team = dash::Team::All();
  size_t nUnits     = team.size();
  
   // cout << "#" << myid << ": " << matOut.pattern().global( {0,0} )[0] << "\n";
   // cout << "#" << myid << ": " << matOut.pattern().local_extents()[0] << endl;
   
   uint gRow =        matOut.pattern().global({0,0})[0];
   uint end  = gRow + matOut.pattern().local_extents()[0];
   double nmax;

   for( uint i = 0; gRow < end; ++gRow, ++i ) {
    nmax = -1;
    for( uint j = 0; j < nelts; ++j ) {
      if( gRow != j) {
        matOut.local[i][j] = distance(points[gRow], points[j]);
        nmax = max( nmax, static_cast<double>( matOut.local[i][j] ) );
      }
    }
    matOut.local[i][gRow] = nelts * nmax;
    vec.local[i] = distance(make_pair(0,0), points[gRow]);
  }
  
  
  dash::barrier( );
  if( 0 == myid ){
    uint count = 0;
    cout << std::showpoint << std::fixed << std::setprecision(5);
    
    for(uint i = 0; i < matOut.extent(0); ++i) {
      for(uint j = 0; j < matOut.extent(1); ++j) {
        if(j) cout << " ";
        cout << static_cast<double>(matOut[i][j]);
      } cout << "\n";
    }
    
    cout << "\n";
    
    for(uint i = 0; i < vec.size(); ++i){
      if(i) cout << " ";
      cout << static_cast<double>(vec[i]);
    } cout << endl;
  }
}
