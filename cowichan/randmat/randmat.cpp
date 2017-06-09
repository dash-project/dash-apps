#include <libdash.h>
#include <iostream>

using std::cout;
using std::endl;
using std::cin;

using uint  = unsigned int ;
using uchar = unsigned char;

struct InputPar { uint nrows, ncols, s;} in;
static int myid;

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

/* 
 * One unit has the job to read in the parameters.
 * Because there's always a unit0, it reads the input parameter and
 * distributes them to the rest of the units.
 */
inline void ReadPars(){

  dash::Shared<InputPar> input_transfer;
  
  if(0 == myid)
  {
    cin >> in.nrows;    
    cin >> in.ncols;
    cin >> in.s;
    
    input_transfer.set(in);
  }
  input_transfer.barrier();
  in = input_transfer.get();
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
void Randmat(T& rand_mat, const uint& nrows, const uint& ncols, const uint& seed)
{
  const int LCG_A = 1664525, LCG_C = 1013904223;

  auto gc   = rand_mat.pattern( ).global( {0,0} );
  uint gbeg = gc[0];  // global row of local (0,0)

  if( 0 < rand_mat.local_size( ) ){
    for( uint i = 0; i < nrows; ++i ) {
      uint s = seed + gbeg + i;

      for( int j = 0; j < ncols; ++j ) {
        s = LCG_A * s + LCG_C;
        rand_mat.lbegin( )[i*ncols + j] = ( (unsigned)s ) % 100;
      }
    }
  }
  dash::barrier( );
}


int main( int argc, char* argv[] )
{
  dash::init( &argc,&argv );

  myid = dash::myid( );
  ReadPars( );

  dash::NArray<unsigned char, 2> rand_mat ( in.nrows, in.ncols );

  Randmat( rand_mat, in.nrows, in.ncols, in.s );
  Print2D( rand_mat );

  dash::finalize( );
}
