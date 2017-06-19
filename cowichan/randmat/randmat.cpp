// #ifndef DASH_ENABLE_LOGGING
  // #define DASH_ENABLE_LOGGING
// #endif

#include <libdash.h>
#include <iostream>

using std::cout;
using std::endl;
using std::cin;

using dash::Shared;

using uint     = unsigned int ;
using uchar    = unsigned char;
using MATRIX_T = uchar        ;

struct InputPar { uint nrows, ncols, s; } in;
static int myid;

#include "randmat.h"

/* 
 * One unit has the job to read in the parameters.
 * Because there's always a unit0, it reads the input parameter and
 * distributes them to the rest of the units.
 */
inline void ReadPars()
{
  Shared<InputPar> input_transfer;
  
  if(0 == myid)
  {
    cin >> in.nrows;    
    cin >> in.ncols;
    cin >> in.s    ;
    
    input_transfer.set(in);
  }
  input_transfer.barrier();
  in = input_transfer.get();
}


/* 
 * This function prints the content of a 2D matrix to std::out.
 * Datatypes are casted to <const uint> for readable output
 * (otherwise uchars would be printed as chars and not as numerics)
 */ 
template< typename T = MATRIX_T >
inline void Print2D( NArray< T, 2 > const & mat )
{
  if(0==myid){
    for( int i = 0; i < mat.extent(0); i++ ) {
      for( int j = 0; j < mat.extent(1); j++ ) {
        cout << std::setw(3) << static_cast<const uint>( mat(i,j) )<< " ";
      }
      cout << endl;
    } cout << endl;
  }
}


int main( int argc, char* argv[] )
{
  dash::init( &argc,&argv );

  myid = dash::myid( );
  ReadPars( );

  NArray<MATRIX_T, 2> rand_mat ( in.nrows, in.ncols );

  Randmat( rand_mat, in.s );
  Print2D( rand_mat );

  dash::finalize( );
}
