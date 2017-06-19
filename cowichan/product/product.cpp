#include <libdash.h>

using std::cin;
using dash::Shared;

using uint  = unsigned int;

uint nelts;
static int myid;

#include "product.h"


inline void ReadMatrixAndVector(
  NArray < double, 2> & matIn,
  vector < double   > & vec  )
{
  if( 0 == myid ) {
    //Read Matrix
    double tmp;
    for ( auto i : matIn ){
      cin >> tmp;
      i = tmp;
    }
    
    //Read Vector
    for (int i = 0; i < vec.size(); i++) {
      cin >> vec[i];
    }
  }
}


inline void ReadNelts( )
{
  Shared<uint> nelts_transfer;

  if(0 == myid)
  {
    cin >> nelts;

    nelts_transfer.set(nelts);
  }
  nelts_transfer.barrier();
  nelts = nelts_transfer.get();
}


int main( int argc, char* argv[] )
{  
  dash::init( &argc,&argv );
  myid = dash::myid( );
  
  ReadNelts( );
  
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