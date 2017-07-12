#include <libdash.h>

using std::cout;
using std::cin;
using std::endl;
using dash::Shared;
using uint  = unsigned int;

uint nelts;
static int myid;

using value = struct{ int row, col; }; //hast to be signed!
#include "outer.h"


inline void PrintOutput(
  NArray < double, 2 > const & matOut ,
  Array  < double    > const & vec    )
{    
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
  
  Shared<uint> nelts_transfer;

  if(0 == myid)
  {
    cin >> nelts;

    nelts_transfer.set(nelts);
  }
  nelts_transfer.barrier();
  nelts = nelts_transfer.get();
}


inline void ReadVectorOfPoints( vector<value> & points ) {
  for( uint i = 0; i < nelts; i++ ) {
    cin >> points[i].row >> points[i].col;
  }
}


int main( int argc, char* argv[] )
{  
  dash::init( &argc,&argv );
  myid = dash::myid( );
  
  ReadNelts( );
  
  vector < value     > points( nelts        );
  NArray < double, 2 > matOut( nelts, nelts );
  Array  < double    > vec   ( nelts        );
  
  //read input points on unit 0 and broadcast to all units
  if( 0 == myid ) ReadVectorOfPoints( points );
  BroadcastPointsToUnits( points );
  
  Outer( points, matOut, vec, nelts );
  PrintOutput(matOut, vec);

  dash::finalize( );
}