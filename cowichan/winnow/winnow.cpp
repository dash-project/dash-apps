#include <libdash.h>

using uint = unsigned int;

// static variables
static struct InputPar { uint nrows, ncols; } in;
static uint   nelts;
static int    myid ;

#include "winnow.h"
#include <time.h>
#include <stdio.h>

std::ifstream raThr_output;

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
    raThr_output.open(argv[1]);
    
    raThr_output >> in.nrows;
    raThr_output >> in.ncols;

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
    int tmp;
    
      for ( auto i : randMat )
      {
        // scanf( "%u", &tmp )  , i = tmp;
        raThr_output >> tmp;
        i = static_cast<T>(tmp);
      }
      
      bool tmpB;
      
      for ( auto i : threshMask )
      {
        // scanf( "%u", &tmpB ) , i = tmpB;
        raThr_output >> tmpB;
        i = static_cast<T>(tmpB);
      }
      
    raThr_output >> nelts;
    raThr_output.close();
    
    nelts_transfer.set(nelts);
  }
  nelts_transfer.barrier();
  nelts = nelts_transfer.get();
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
  ReadRowsNCols( argv );
  
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
  
  if( clock_gettime( CLOCK_MONOTONIC_RAW, &start) == -1 ) {
    perror( "clock gettime error 1" );
    exit( EXIT_FAILURE );
  }
  
  Winnow( in.nrows, in.ncols, randMat, threshMask, nelts, result );
  
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
    fprintf( fp, "DASH,Winnow,%u, %u, , %u, %u, %.9lf\n", in.nrows, in.ncols, nelts, dash::Team::All().size(), accum );
    fclose ( fp );
  }
  
  // output
  if( 0 == myid )
  {
    cout << nelts << "\n";
    
    for( value it : result ) cout << it;
    
    cout << endl;
  }
  dash::finalize( );
}


