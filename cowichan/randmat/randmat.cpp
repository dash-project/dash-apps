#include <libdash.h>
#include <iostream>

using std::cout;
using std::endl;
using std::cin;

using uint  = unsigned int ;
using uchar = unsigned char;

struct inputPar { uint nrows, ncols, s;};

/* 
 * This function prints the content of a 2D matrix to std::out.
 * Datatypes are casted to <const uint> for readable output
 * (otherwise uchars would be printed as chars and not as numerics)
 */ 
template< typename T >
inline void print2d(const T& mat ) {
  for( int i = 0; i < mat.extent(0); i++ ) {
    for( int j = 0; j < mat.extent(1); j++ ) {
      cout << std::setw(3) << static_cast<const uint>( mat(i,j) )<< " ";
    }
    cout << endl;
  }
}

/* 
 * One unit has the job to read in the parameters.
 * Because there's always a unit0, it reads the input parameter and
 * distributes them to the rest of the units.
 */
template<typename T>
inline void readPars(const T myid, inputPar& input){

  dash::Shared<inputPar> inputTransfer;
  
  if(0 == myid)
  {
    cin >> input.nrows;    
    cin >> input.ncols;
    cin >> input.s;
    
    inputTransfer.set(input);
  }
  inputTransfer.barrier();
  input = inputTransfer.get();
}

/*
 * The Cowichan problems require that the output is independent of the
 * numbers of processors used. For randmat() a common solution found
 * in other implementations is to use a simple linear congruential
 * random number generator (LCG) with a separate deterministic seed
 * for each row and to parallelize over the rows of the matrix. This
 * is also how the DASH solution below works.
 *
 * A potential alternative would be to use a counter-based random
 * number generation scheme (e.g. random123) that can be easily
 * parallelized.
 */

template< typename T >
void randmat( dash::NArray<T, 2>& mat, const uint& nrows, const uint& ncols, const uint& seed )
{
  const int LCG_A = 1664525, LCG_C = 1013904223;

  auto gc   = mat.pattern( ).global( {0,0} );
  uint gbeg = gc[0];  // global row of local (0,0)

  if( 0 < mat.local_size( ) ){
    for( uint i = 0; i < nrows; ++i ) {
      uint s = seed + gbeg + i;

      for( int j = 0; j < ncols; ++j ) {
        s = LCG_A * s + LCG_C;
        mat.lbegin( )[i*ncols + j] = ( (unsigned)s ) % 100;
      }
    }
  }
}


int main( int argc, char* argv[] )
{
  dash::init( &argc,&argv );
  auto myid = dash::myid( );

  inputPar input;
  readPars(myid, input);

  dash::NArray<unsigned char, 2> rand_mat ( input.nrows, input.ncols );

  randmat( rand_mat, input.nrows, input.ncols, input.s );

  dash::barrier( );
  if( myid==0 ) print2d( rand_mat );

  dash::finalize( );
}
