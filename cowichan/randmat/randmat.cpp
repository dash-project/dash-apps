#include <libdash.h>
#include <iostream>

using std::cout;
using std::endl;
using std::cin;

using uint  = unsigned int ;
using uchar = unsigned char;


template< typename T >
inline void print2d(const T& mat ) {
  for( int i = 0; i < mat.extent(0); i++ ) {
    for( int j = 0; j < mat.extent(1); j++ ) {
      cout << std::setw(3) << static_cast<const uint>( mat(i,j) )<< " ";
    }
    cout << endl;
  }
}

inline void readPars( uint &nrows, uint &ncols, uint &s){
  dash::Shared<uint> nrows, ncols, s;

  uint tmp;
  cin >> tmp;
  nrows.set(tmp);
  cin >> tmp;
  ncols.set(tmp);
  cin >> tmp;
  s.set(tmp);
  
  nrows.flush();
  ncols.flush();
  s.flush();
  }
  dash::barrier();
}

//
// The Cowichan problems require that the output is independent of the
// numbers of processors used. For randmat() a common solution found
// in other implementations is to use a simple linear congruential
// random number generator (LCG) with a separate deterministic seed
// for each row and to parallelize over the rows of the matrix. This
// is also how the DASH solution below works.
//
// A potential alternative would be to use a counter-based random
// number generation scheme (e.g. random123) that can be easily
// parallelized.
//

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

  uint nrows, ncols, s;
  readPars( nrows, ncols, s);
  
  cout << "test\n" << "nrows:" << nrows << ", ncols:" << ncols << endl;
  

  dash::NArray<unsigned char, 2> rand_mat ( nrows, ncols );

  randmat( mat, nrows, ncols, s );

  dash::barrier( );
  if( myid==0 ) print2d( rand_mat );

  dash::finalize( );
}
