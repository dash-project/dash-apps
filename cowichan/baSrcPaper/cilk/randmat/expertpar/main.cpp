/*
 * randmat: random number generation
 *
 * input:
 *   nrows, ncols: the number of rows and columns
 *   s: the seed
 *
 * output:
 *   matrix: an nrows x ncols integer matrix
 */

#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
//export CILK_NWORKERS=24 && echo "50000 42950 55" | ./main --is_bench --nproc 24
static int *matrix;

int is_bench = 0;

void randmat(int nrows, int ncols, int seed) {
  const int LCG_A = 1664525, LCG_C = 1013904223;
  cilk_for (int i = 0; i < nrows; i++) {
    int begin = i;
    int s = seed + i;
    for (int j = 0; j < ncols; j++) {
      s = LCG_A * s + LCG_C;
      matrix[begin*ncols + j] = ((unsigned)s) % 100;
    }
  }
}

int main(int argc, char *argv[]) {
  int nrows, ncols, s, i, j, nproc = 0;
  struct timespec start, stop;
  double accum;

  if (argc >= 2) {
    for (int a = 0; a < argc; a++) {
      if (!strcmp(argv[a], "--is_bench")) {
        is_bench = 1;
      } else if (!strcmp(argv[a], "--nproc")) {
        sscanf(argv[++a], "%d", &nproc);
      }
    }
  }

  scanf("%d%d%d", &nrows, &ncols, &s);
  matrix = (int*) malloc (sizeof(int) * nrows * ncols);
  
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
  fprintf( fp, "Cilk  ,Randmat,%5u,%5u,   ,     ,%2u,%.9lf,isBench:%d\n", nrows, ncols, nproc, accum, is_bench );
  fclose ( fp );

  
  if (!is_bench) {
    for (i = 0; i < nrows; i++) {
      for (j = 0; j < ncols; j++) {
        printf("%d ", matrix[i*ncols + j]);
      }
      printf("\n");
    }
    printf("\n");
  }
  
  return 0;
}
