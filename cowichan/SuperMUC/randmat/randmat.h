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

using dash::barrier;
using dash::NArray;

template< typename T = MATRIX_T >
inline void Randmat(
  NArray< T, 2 >       & rand_mat,
            uint const   seed    )
{
  const uint LCG_A = 1664525, LCG_C = 1013904223;

  uint nrows = rand_mat.local.extent(0); // num of local rows
  uint ncols = rand_mat.local.extent(1); // num of local cols

  auto gc   = rand_mat.pattern( ).global( {0,0} );
  uint gbeg = gc[0];  // global row of local (0,0)

  if( 0 < rand_mat.local_size( ) ){
    for( uint i = 0; i < nrows; ++i ) {
      uint s = seed + gbeg + i;

      for( uint j = 0; j < ncols; ++j ) {
        s = LCG_A * s + LCG_C;
        rand_mat.lbegin( )[i*ncols + j] = ( (unsigned)s ) % 1000;
      }
    }
  }
  barrier( );
}