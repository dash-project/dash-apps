#include <libdash.h>
#include <iostream>
#include <cstdio>

using std::cout;
using std::endl;
using std::cin;
using uint     = unsigned int ;
using uchar    = unsigned char;
using MATRIX_T = uchar        ;

struct     InputPar { uint nrows, ncols; } in;
uint       percent;
static int myid;

#include "thresh.h"
#include <time.h>
#include <stdio.h>

std::ifstream randmat_output;

/*
 * This function prints the content of a 2D matrix to std::out.
 * Datatypes are casted to <const uint> for readable output
 * (otherwise uchars would be printed as chars and not as numerics)
 */
inline void Print2D( NArray< bool, 2 > const & mat )
{
  if(0==myid){

    //cout << in.nrows << " " << in.ncols << "\n";

    for( int i = 0; i < mat.extent(0); i++ ) {
      for( int j = 0; j < mat.extent(1); j++ ) {
        cout << static_cast<const uint>( mat(i,j) )<< " ";
      }
      cout << "\n";
    }
  }
}

// unit0 reads the matrix with random values from std::in
template<typename T = MATRIX_T>
inline void ReadRandMat( NArray< T, 2 > & rand_mat )
{
  if(0 == myid){
    int tmp;
    for ( auto i : rand_mat ){
      randmat_output >> tmp;
      i = static_cast<T>(tmp);
    }
  }
}

/*
 * One unit has the job to read in the parameters.
 * Because there's always a unit0, it reads the input parameter and
 * distributes them to the rest of the units.
 */
inline void ReadRowsNCols( char * argv[] )
{
  Shared<InputPar> input_transfer;

  if(0 == myid)
  {
    randmat_output.open(argv[1]);
    
    randmat_output >> in.nrows;
    randmat_output >> in.ncols;

    input_transfer.set(in);
  }
  input_transfer.barrier();
  in = input_transfer.get();
}

inline void ReadPercentage( )
{
  Shared<uint> percent_transfer;

  if(0 == myid)
  {
    randmat_output >> percent;

    percent_transfer.set(percent);
    
    randmat_output.close();
  }
  percent_transfer.barrier();
  percent = percent_transfer.get();
}


int main( int argc, char* argv[] )
{
  dash::init( &argc, &argv );

  struct timespec start, stop;
  double accum;
  int is_bench = 0;
  
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--is_bench")) {
      is_bench = 1;
    }
  }
  
  myid = dash::myid( );
  ReadRowsNCols( argv );

  NArray< MATRIX_T, 2 > rand_mat    ( in.nrows, in.ncols );
  NArray< bool    , 2 > thresh_mask ( in.nrows, in.ncols );

  ReadRandMat(rand_mat);
  ReadPercentage();
  
  if( clock_gettime( CLOCK_MONOTONIC_RAW, &start) == -1 ) {
    perror( "clock gettime error 1" );
    exit( EXIT_FAILURE );
  }
  
  Thresh( rand_mat, thresh_mask, in.nrows, in.ncols, percent );
  
  if( clock_gettime( CLOCK_MONOTONIC_RAW, &stop) == -1 ) {
    perror( "clock gettime error 2" );
    exit( EXIT_FAILURE );
  }
  
  accum = ( stop.tv_sec - start.tv_sec ) + ( stop.tv_nsec - start.tv_nsec ) / 1e9;
  

  if( is_bench && 0 == myid ){
    FILE* fp = fopen("./measurements.txt", "a");
    
    if( !fp ) {
        perror("File opening for benchmark results failed");
        return EXIT_FAILURE;
    }
    // Lang, Problem, rows, cols, thresh, winnow_nelts, jobs, time
    fprintf( fp, "DASH,Thresh,%u, %u, %u , , %u, %.9lf\n", in.nrows, in.ncols, percent, dash::Team::All().size(), accum );
    fclose ( fp );
  }
  
  Print2D( thresh_mask );
  dash::finalize( );
}
