#include <libdash.h>

using std::cin;
using dash::Shared;

using uint  = unsigned int;

uint nelts;
static int myid;

#include "product.h"

std::ifstream outer_output;

inline void ReadMatrixAndVector(
  NArray < double, 2> & matIn,
  vector < double   > & vec  )
{
  if( 0 == myid )
  {
    // read matrix
    double tmp;
    for ( auto i : matIn ){ outer_output >> tmp; i = tmp; }
    
    // //Read Vector
    for (int i = 0; i < vec.size(); i++){ outer_output >> vec[i]; }
    
    outer_output.close();
  }
}


inline void ReadNelts( char * argv[] )
{
  Shared<uint> nelts_transfer;

  if(0 == myid)
  {
    outer_output.open(argv[1]);
    outer_output >> nelts;

    nelts_transfer.set(nelts);
  }
  nelts_transfer.barrier();
  nelts = nelts_transfer.get();
}


int main( int argc, char* argv[] )
{  
  dash::init( &argc,&argv );
  myid = dash::myid( );
  
  ReadNelts( argv );

  NArray < double, 2 > matIn  ( nelts, nelts );
  Array  < double    > result ( nelts        );
  vector < double    > vec    ( nelts        );
 
  //read input on unit 0 and broadcast it
  ReadMatrixAndVector(matIn, vec);
  BroadcastOuterVecToUnits(vec);

  Product(vec, matIn, result, nelts );
  PrintOutput( result, nelts );
  
  dash::finalize( );
}