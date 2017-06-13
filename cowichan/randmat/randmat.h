/* These "usings" and definitionswill be also in chain! */
using std::cout;
using std::endl;
using std::cin;

using uint  = unsigned int ;
using uchar = unsigned char;

struct InputPar { uint nrows, ncols, s;} in;
static int myid;

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
inline void Randmat(T& rand_mat, const uint nrows, const uint ncols, const uint seed)
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