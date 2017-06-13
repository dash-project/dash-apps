/* #ifndef DASH_ENABLE_LOGGING
#define DASH_ENABLE_LOGGING
#endif */

#include <libdash.h>

using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::pair;
using std::max;
using std::make_pair;

using uint  = unsigned int;
using POI_T = int;  //this type musst be signed!


uint nelts;
static int myid;
#include "outer.h"


template< typename T, typename Y >
inline void PrintOutput(T const& matOut, Y const& vec ) {
  if( 0 == myid ){
    cout << nelts << "\n";
    uint count = 0;
    cout << std::showpoint << std::fixed << std::setprecision(4);
    
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


inline void ReadNelts( ){
  
  dash::Shared<uint> nelts_transfer;

  if(0 == myid)
  {
    cin >> nelts;

    nelts_transfer.set(nelts);
  }
  nelts_transfer.barrier();
  nelts = nelts_transfer.get();
}


template<typename T = POI_T>
inline void ReadVectorOfPoints( vector< pair<T,T> > & points ) {
  for( uint i = 0; i < nelts; i++ ) {
    cin >> points[i].first >> points[i].second;
  }
}


template<typename T = POI_T>
inline void BroadcastInputToUnits( vector< pair<T,T> > & points ) {
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
  myid = dash::myid( );
  
  ReadNelts( );
  
  vector< pair<POI_T, POI_T> > points( nelts );
  dash::NArray< double, 2 >    matOut( nelts, nelts );
  dash::Array < double >       vec( nelts );
  
  //read input points on unit 0 and broadcast to all units
  if( 0 == myid ) ReadVectorOfPoints( points );
  BroadcastInputToUnits( points );
  
  Outer( points, matOut, vec );
  PrintOutput(matOut, vec);

  dash::finalize( );
}