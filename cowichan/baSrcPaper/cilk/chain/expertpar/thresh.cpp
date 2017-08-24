/*
 * thresh: histogram thresholding
 *
 * input:
 *   randmat_matrix: the integer matrix to be thresholded
 *   nrows, ncols: the number of rows and columns
 *   percent: the percentage of cells to retain
 *
 * output:
 *   thresh_mask: a boolean matrix filled with true for cells kept
 */

#include <cilk/cilk.h>
#include <cilk/reducer_max.h>
#include <cilk/cilk_api.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <algorithm>

extern int *randmat_matrix;
int *thresh_mask;
static int histogram[128][200];


int reduce_max (int nrows, int ncols) {
  cilk::reducer_max <int> max_reducer (0);

  cilk_for (int i = 0; i < nrows; i++) {
    int begin = i;
    int tmp_max = 0;
    for (int j = 0; j < ncols; j++) {
      tmp_max = std::max (tmp_max, randmat_matrix [begin*ncols + j]);
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
      histogram [Self][randmat_matrix[r*ncols +i]]++;
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

void fill_thresh_mask (int nrows, int ncols, int threshold) {
  cilk_for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
      thresh_mask[i*ncols +j] = randmat_matrix [i*ncols +j] >= threshold;
    }
  }
}

void thresh(int nrows, int ncols, int percent) {
  int i;
  int nmax = 0;
  int count, prefixsum, threshold;

  thresh_mask = (int*) malloc (sizeof(int) * ncols * nrows);

  nmax = reduce_max(nrows, ncols);

  fill_histogram(nrows, ncols);
  merge_histogram();

  count = (static_cast<size_t>(nrows) * ncols * percent) / 100;
  // printf("threshCount%i\n", count);

  prefixsum = 0;
  threshold = nmax;

  for (i = nmax; i >= 0 && prefixsum <= count; i--) {
    prefixsum += histogram[0][i];
    threshold = i;
  }

  fill_thresh_mask(nrows, ncols, threshold);
}
