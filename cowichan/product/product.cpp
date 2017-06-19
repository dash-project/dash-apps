/* #ifndef DASH_ENABLE_LOGGING
#define DASH_ENABLE_LOGGING
#endif */

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


inline void BroadcastInputToUnits( vector <double> & vec )
{
  team_unit_t TeamUnit0ID = Team::All().myid( );
  TeamUnit0ID.id = 0;
  dart_ret_t ret = dart_bcast(
                      static_cast<void*>(vec.data( )),  // buf 
                      vec.size( )                    ,  // nelts
                      DART_TYPE_DOUBLE               ,  // dtype
                      TeamUnit0ID                    ,  // root
                      Team::All().dart_id( )         ); // team
                      
  if( DART_OK != ret ) cout << "An error while BCAST has occured!" << endl; 
}


int main( int argc, char* argv[] )
{  
  dash::init( &argc,&argv );
  myid = dash::myid( );
  
  ReadNelts( );
  
  NArray < double, 2 > matIn  ( nelts, nelts);
  Array  < double    > result ( nelts       );
  vector < double    > vec    ( nelts       );
 
  //read input on unit 0 and broadcast it
  ReadMatrixAndVector(matIn, vec);
  BroadcastInputToUnits(vec);

  Product(vec, matIn, result);
  PrintOutput( result );
  
  dash::finalize( );
}