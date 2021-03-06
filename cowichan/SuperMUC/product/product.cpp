#include <libdash.h>

using std::cin;
using dash::Shared;

using uint  = unsigned int;

uint nelts;
static int myid;

#include "product.h"
#include <cstring>
#include <time.h>
#include <stdio.h>

using std::strcmp;

std::ifstream outer_output;

inline void ReadMatrixAndVector(
  NArray < double, 2> & matIn,
  vector < double   > & vec  )
{
  if( 0 == myid )
  {
    // read matrix
    double tmp;
    for ( auto i : matIn ){ outer_output >> tmp; i = tmp; }

    // //Read Vector
    for (int i = 0; i < vec.size(); i++){ outer_output >> vec[i]; }

    //outer_output.close();
  }
}


inline void ReadNelts( char * argv[] )
{
  Shared<uint> nelts_transfer;

  if(0 == myid)
  {
    //changed to cin again for SuperMUC benchmark
    //outer_output.open(argv[1]);
    cin >> nelts;

    nelts_transfer.set(nelts);
  }
  nelts_transfer.barrier();
  nelts = nelts_transfer.get();
}


int main( int argc, char* argv[] )
{
  dash::init( &argc,&argv );
  myid = dash::myid( );

  struct timespec start, stop;
  double accum;
  int is_bench = 0;

  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--is_bench")) {
      is_bench = 1;
    }
  }

  ReadNelts( argv );

  NArray < double, 2 > matIn  ( nelts, nelts );
  Array  < double    > result ( nelts        );
  vector < double    > vec    ( nelts        );

  //read input on unit 0
  if (!is_bench) { ReadMatrixAndVector(matIn, vec); }

  if( clock_gettime( CLOCK_MONOTONIC, &start) == -1 ) {
    perror( "clock gettime error 1" );
    exit( EXIT_FAILURE );
  }

  //broadcast vector from unit0 to all other units
  BroadcastOuterVecToUnits(vec);
  Product(vec, matIn, result, nelts );

  if( clock_gettime( CLOCK_MONOTONIC, &stop) == -1 ) {
    perror( "clock gettime error 2" );
    exit( EXIT_FAILURE );
  }

  accum = ( stop.tv_sec - start.tv_sec ) + ( stop.tv_nsec - start.tv_nsec ) / 1e9;


  if( 0 == myid ){
    // Lang, Problem, rows, cols, thresh, winnow_nelts, jobs, time
    printf( "DASH  ,Product,     ,     ,   ,%5u,%2u,%.9lf,isBench:%d\n", nelts, dash::Team::All().size(), accum, is_bench );
  }

  if( !is_bench ){ PrintOutput( result, nelts ); }
  dash::finalize( );
}