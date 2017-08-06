/*
 * chain: chain all problems
 *
 * input:
 *   nelts: the number of elements
 *   randmat_seed: random number generator of cells to retain
 *   thresh_percent: percentage of cells to retain
 *   winnow_nelts: the number of points to select
 *
 * output:
 *   result: a real vector, whose values are the result of the final product
 */
#include <cilk/cilk.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "common.h"

int is_bench = 0;

extern double *product_result;

int main(int argc, char** argv) {
  int nelts, randmat_seed, thresh_percent, winnow_nelts, i, nproc = 0;
  struct timespec start, stop;
  double accum;

  if (argc >= 2) {
    for (int a = 0; a < argc; a++){
      if (!strcmp(argv[a], "--is_bench")) {
        is_bench = 1;
      } else if (!strcmp(argv[a], "--nproc")) {
        sscanf(argv[++a], "%d", &nproc);
      }
    }
  }

  scanf("%d%d%d%d", &nelts, &randmat_seed, &thresh_percent, &winnow_nelts);
  
  if( clock_gettime( CLOCK_MONOTONIC_RAW, &start) == -1 ) {
    perror( "clock gettime error 1" );
    exit( EXIT_FAILURE );
  }
  
  randmat(nelts, nelts, randmat_seed); cilk_sync;
  thresh(nelts, nelts, thresh_percent); cilk_sync;
  winnow(nelts, nelts, winnow_nelts); cilk_sync;
  outer(winnow_nelts); cilk_sync;
  product(winnow_nelts); cilk_sync;
  
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
  fprintf( fp, "Cilk,Chain, %u, %u, %u, %u, %u, %.9lf,isBench:%d\n", nelts, nelts, thresh_percent, winnow_nelts, nproc, accum, is_bench );
  fclose ( fp );
  

  if (!is_bench) {
    printf("%d\n", winnow_nelts);
    for (i = 0; i < winnow_nelts; i++) {
      printf("%.4f ", product_result[i]);
    }
    printf("\n");
  }

  return 0;
}

