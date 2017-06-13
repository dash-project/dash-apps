#include <libdash.h>
#include <iostream>

#define MATRIX_T uchar
#include "randmat.h"

/* 
 * This function prints the content of a 2D matrix to std::out.
 * Datatypes are casted to <const uint> for readable output
 * (otherwise uchars would be printed as chars and not as numerics)
 */ 
template< typename T >
inline void Print2D(const T& mat ) {
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

  dash::NArray<MATRIX_T, 2> rand_mat ( in.nrows, in.ncols );

  Randmat( rand_mat, in.nrows, in.ncols, in.s );
  Print2D( rand_mat );

  dash::finalize( );
}
