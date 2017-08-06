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
#include <time.h>
#include <stdio.h>
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

  struct timespec start, stop;
  double accum;
  int is_bench = 0;
  
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--is_bench")) {
      is_bench = 1;
    }
  }
  
  myid = dash::myid( );
  ReadPars( );

  NArray<MATRIX_T, 2> rand_mat ( in.nrows, in.ncols );

  if( clock_gettime( CLOCK_MONOTONIC_RAW, &start) == -1 ) {
    perror( "clock gettime error 1" );
    exit( EXIT_FAILURE );
  }
  
  Randmat( rand_mat, in.s );
  
  if( clock_gettime( CLOCK_MONOTONIC_RAW, &stop) == -1 ) {
    perror( "clock gettime error 2" );
    exit( EXIT_FAILURE );
  }
  
  accum = ( stop.tv_sec - start.tv_sec ) + ( stop.tv_nsec - start.tv_nsec ) / 1e9;
  
  
  if( 0 == myid ){
    FILE* fp = fopen("./measurements.txt", "a");
    
    if( !fp ) {
        perror("File opening for benchmark results failed");
        return EXIT_FAILURE;
    }
    // Lang, Problem, rows, cols, thresh, winnow_nelts, jobs, time
    fprintf( fp, "DASH,Randmat,%u, %u, , , %u, %.9lf,isBench:%d\n", in.nrows, in.ncols, dash::Team::All().size(), accum, is_bench );
    fclose ( fp );
  }

  if (!is_bench) { Print2D( rand_mat ); }
  dash::finalize( );
}