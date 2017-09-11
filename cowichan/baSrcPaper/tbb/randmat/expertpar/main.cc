/* randmat: random number generator
 *
 * input:
 *   nrows, ncols: number of rows and columns
 *   s: random number generation seed
 *
 * output:
 *   matrix: random nrows x ncols integer matrix
 */
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/task_scheduler_init.h"
#include <time.h>
#include <stdio.h>

using namespace tbb;

typedef blocked_range<size_t> range;

static unsigned char* matrix;

int is_bench = 0;
int n_threads = task_scheduler_init::default_num_threads();

void randmat(int nrows, int ncols, unsigned int s) {
  const int LCG_A = 1664525, LCG_C = 1013904223;
  parallel_for(
    range(0, nrows),
    [=](range r) {
      auto end = r.end (); 
      for (size_t i = r.begin(); i != end; ++i) {
        unsigned int seed = s + i;
        for (int j = 0; j < ncols; j++) {
          seed = LCG_A * seed + LCG_C;
          matrix[i*ncols + j] = seed % 100;
        }
      }
  });
}

int main(int argc, char** argv) {
  int nrows, ncols, s;
  struct timespec start, stop;
  double accum;

  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--is_bench")) {
      is_bench = 1;
    } else if (!strcmp(argv[i], "--threads")) {
      sscanf(argv[i + 1], "%d", &n_threads);
      i++;
    }
  }

  task_scheduler_init init(n_threads);

  scanf("%d%d%d", &nrows, &ncols, &s);
  matrix = (unsigned char*) malloc (sizeof (unsigned char) * nrows * ncols);
  
  if( clock_gettime( CLOCK_MONOTONIC_RAW, &start) == -1 ) {
    perror( "clock gettime error 1" );
    exit( EXIT_FAILURE );
  }
  
  randmat(nrows, ncols, s);
  
  if( clock_gettime( CLOCK_MONOTONIC_RAW, &stop) == -1 ) {
    perror( "clock gettime error 2" );
    exit( EXIT_FAILURE );
  }
  
  accum = ( stop.tv_sec - start.tv_sec ) + ( stop.tv_nsec - start.tv_nsec ) / 1e9;
  
  

  FILE* fp = fopen("./measurements.txt", "a");
  
  if( !fp ) {
      perror("File opening for benchmark results failed");
      return EXIT_FAILURE;
  }
  // Lang, Problem, rows, cols, thresh, winnow_nelts, jobs, time
  fprintf( fp, "TBB   ,Randmat,%5u,%5u,   ,     ,%2u,%.9lf, isBench:%d\n", nrows, ncols, n_threads, accum, is_bench );
  fclose ( fp );
  

  if (!is_bench) {
    //printf("%d %d\n", nrows, ncols);
    for (int i = 0; i < nrows; i++) {
      for (int j = 0; j < ncols; j++) {
        printf("%d ", matrix[i*ncols + j]);
      }
      printf("\n");
    }
    printf("\n");
  }

  return 0;
}
