#include <libdash.h>

using uint = unsigned int;

// static variables
static struct InputPar { uint nrows, ncols; } in;
static uint   nelts;
static int    myid ;

#include "winnow.h"


/*
 * One unit has the job to read in the parameters.
 * Because there's always a unit0, it reads the input parameter and
 * distributes them to the rest of the units.
 */
inline void ReadRowsNCols( )
{
  Shared<InputPar> input_transfer;

  if(0 == myid)
  {
    cin >> in.nrows;
    cin >> in.ncols;

    input_transfer.set(in);
  }
  input_transfer.barrier();
  in = input_transfer.get();
}
 

template< typename T = MATRIX_T >
inline void ReadMatricesAndNelts( NArray<T,2>& randMat, NArray<bool,2>& threshMask )
{
  Shared<uint> nelts_transfer;

  if(0 == myid)
  {
    //read matrices
    T tmp;
      for ( auto i : randMat ){
        scanf( "%u", &tmp )  , i = tmp;
      }
      bool tmpB;
      for ( auto i : threshMask ){
        scanf( "%u", &tmpB ) , i = tmpB;
      }
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
  
  ReadRowsNCols( );
  
  NArray< MATRIX_T, 2 > randMat    ( in.nrows, in.ncols );
  NArray< bool    , 2 > threshMask ( in.nrows, in.ncols );
  
  res_array_t result;
  

  #ifdef DEBUG  // print error message if mask's and matrix's local size aren't identical
    if( threshMask.local_size() != randMat.local_size() )
    {
      cout << "On unit " << myid
           << " the local sizes of matrix and mask differ!\naborted on this unit\n";
      return -1;
    }
  #endif
  
  ReadMatricesAndNelts( randMat, threshMask );
  
  
  Winnow( in.nrows, in.ncols, randMat, threshMask, nelts, result );
  
  // output
  if( 0 == myid )
  {
    cout << nelts << "\n";
    
    for( value it : result ) cout << it;
    
    cout << endl;
  }
  
  dash::finalize( );
}


