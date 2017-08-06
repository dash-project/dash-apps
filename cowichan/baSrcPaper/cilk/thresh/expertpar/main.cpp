/*
 * thresh: histogram thresholding
 *
 * input:
 *   matrix: the integer matrix to be thresholded
 *   nrows, ncols: the number of rows and columns
 *   percent: the percentage of cells to retain
 *
 * output:
 *   mask: a boolean matrix filled with true for cells kept
 */

#include <cilk/cilk.h>
#include <cilk/reducer_max.h>
#include <cilk/cilk_api.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <algorithm>

int is_bench = 0;
static int *matrix;
static int *mask;

static int histogram[128][200];

int reduce_max (int nrows, int ncols) {
  cilk::reducer_max <int> max_reducer (0);

  cilk_for (int i = 0; i < nrows; i++) {
    int begin = i;
    int tmp_max = 0;
    for (int j = 0; j < ncols; j++) {
      tmp_max = std::max (tmp_max, matrix [begin*ncols + j]);
    }
    max_reducer.calc_max (tmp_max);
  }

  return max_reducer.get_value ();
}

void fill_histogram(int nrows, int ncols) {
  int P = __cilkrts_get_nworkers();
  cilk_for (int r = 0; r < nrows; ++r) {
    int Self = __cilkrts_get_worker_number();
    for (int i = 0; i < ncols; i++) {
      histogram [Self][matrix[r*ncols + i]]++;
    }
  }
}

void merge_histogram () {
  int P = __cilkrts_get_nworkers();
  cilk_for (int v = 0; v < 100; ++v) {
    int merge_val = __sec_reduce_add (histogram [1:(P-1)][v]);
    histogram [0][v] += merge_val;
  }
}

void fill_mask (int nrows, int ncols, int threshold) {
  cilk_for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
      mask[i*ncols + j] = matrix [i*ncols + j] >= threshold;
    }
  }
}

void thresh(int nrows, int ncols, int percent) {
  int i;
  int nmax = 0;
  int count, prefixsum, threshold;

  nmax = reduce_max(nrows, ncols);

  fill_histogram(nrows, ncols);
  merge_histogram();

  count = (nrows * ncols * percent) / 100;

  prefixsum = 0;
  threshold = nmax;

  for (i = nmax; i >= 0 && prefixsum <= count; i--) {
    prefixsum += histogram[0][i];
    threshold = i;
  }

  fill_mask(nrows, ncols, threshold);
}

int main(int argc, char *argv[]) {
  int nrows, ncols, percent, i, j, nproc = 0;
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

  scanf("%d%d", &nrows, &ncols);
  matrix = (int*) malloc (sizeof(int) * ncols * nrows);
  mask = (int*) malloc (sizeof(int) * ncols * nrows);

  if (!is_bench) {
    for (i = 0; i < nrows; i++) {
      for (j = 0; j < ncols; j++) {
        scanf("%hhu", &matrix[i*ncols + j]);
      }
    }
  }

  scanf("%d", &percent);
  
  if( clock_gettime( CLOCK_MONOTONIC_RAW, &start) == -1 ) {
    perror( "clock gettime error 1" );
    exit( EXIT_FAILURE );
  }

  thresh(nrows, ncols, percent);
  
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
  fprintf( fp, "Cilk,Thresh,%u, %u, %u, , %u, %.9lf,isBench:%d\n", nrows, ncols, percent, nproc, accum, is_bench );
  fclose ( fp );
  

  if (!is_bench) {
    //printf("%d %d\n", nrows, ncols);
    for (i = 0; i < nrows; i++) {
      for (j = 0; j < ncols; j++) {
        printf("%hhu ", mask[i*ncols + j]);
      }
      printf("\n");
    }
    printf("\n");
  }

  return 0;
}
