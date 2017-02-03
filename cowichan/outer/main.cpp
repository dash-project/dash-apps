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
using std::this_thread::sleep_for;
using std::pair;

using uint       = unsigned  int;
using POI_T = uint;

template<typename T = POI_T>
inline void outer( uint nelts, vector< pair<T,T> > & points, dash::NArray< double, 2 > & matOut, int myid );

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
                      points.size( ) * sizeof(pair<T,T>)  ,  // nelem
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
  
  vector< pair<POI_T, POI_T> > points(nelts);
  dash::NArray< double, 2 >    matOut( nelts, nelts );
  
  //read input points on unit 0
  if( 0 == myid ) read_vector_of_points( nelts, points );
  
  //broadcast input from unit 0
  broadcastInputToUnits( nelts, points );

  outer( nelts, points, matOut, myid );

  dash::finalize( );
}


template<typename T = POI_T>
inline void outer( uint nelts, vector< pair<T,T> > & points, dash::NArray< double, 2 > & matOut, int myid ){

  dash::Team & team = dash::Team::All();
  size_t nUnits     = team.size();

  
  
}
