#include <libdash.h>

using uint = unsigned int;

// static variables
static struct InputPar { uint nrows, ncols; } in;
static uint   nelts;
static uint   thresh;
static int    myid ;

#include "winnow.h"
#include <time.h>
#include <stdio.h>

std::ifstream raThr_output;
  int is_bench = 0;

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
  Shared<uint> thresh_transfer;

  if(0 == myid)
  {
    if (!is_bench) { 
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
    }
      
    raThr_output >> nelts;
    if(is_bench) raThr_output >> thresh;
    raThr_output.close();
    
    thresh_transfer.set(thresh);
    nelts_transfer.set(nelts);
  }
  nelts_transfer.barrier();
  nelts = nelts_transfer.get();
  thresh = thresh_transfer.get();
}

template< typename T = MATRIX_T >
inline void FillOnTheFly( NArray<T,2>& randMat, NArray<bool,2>& threshMask )
{
  auto gR = randMat.pattern( ).global( {0,0} );
  auto gT = threshMask.pattern( ).global( {0,0} );
  
  uint i = gR[0];  // global row of local (0,0)
  uint j = 0;
  // cout << "ich hab:" << i << "\n";
  
  for( T * ptr = randMat.lbegin(); ptr < randMat.lend(); ++ptr)
  {
    *ptr = (i*in.ncols+j) % 100;
    if(++j == in.ncols) {++i;j=0;}
  }
  
  uint threshInverse = 100 / thresh;
  i =  gT[0];  // global row of local (0,0)
  j = 0;
  
  for( bool * ptr = threshMask.lbegin(); ptr < threshMask.lend(); ++ptr) 
  {
    if(( ( i*in.ncols+j) % threshInverse ) == 0) {
      *ptr = true;
    }else{
      *ptr = false;
    }
    if(++j == in.ncols) {++i;j=0;}
  }
   
  dash::barrier();
}


int main( int argc, char* argv[] )
{  
  dash::init( &argc,&argv );
  
  struct timespec start, stop;
  double accum;
  
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
  
  if (is_bench) FillOnTheFly( randMat, threshMask );
  
  if( clock_gettime( CLOCK_MONOTONIC, &start) == -1 ) {
    perror( "clock gettime error 1" );
    exit( EXIT_FAILURE );
  }
  
  Winnow( in.nrows, in.ncols, randMat, threshMask, nelts, result );
  
  if( clock_gettime( CLOCK_MONOTONIC, &stop) == -1 ) {
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
    fprintf( fp, "DASH  ,Winnow ,%5u,%5u,%3u,%5u,%2u,%.9lf,isBench:%d\n", in.nrows, in.ncols, thresh, nelts, dash::Team::All().size(), accum, is_bench );
    fclose ( fp );
  }
  
  if (!is_bench) { 
    // output
    if( 0 == myid )
    {
      cout << nelts << "\n";
      
      for( value it : result ) cout << it;
      
      cout << endl;
    }
  }
  dash::finalize( );
}


