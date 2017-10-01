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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

extern int *randmat_matrix;
int *thresh_mask;
static int histogram[128][100];

int max(int a, int b) {
  return a > b ? a : b;
}

// parallel max reduce on [begin, end)
int reduce_max(int begin, int end, int ncols) {
  int middle = begin + (end - begin) / 2;
  int left, right, res, i;
  if (begin + 1 == end) {
    res = randmat_matrix[begin*ncols +0];
    for (i = 1; i < ncols; i++) {
      res = max(res, randmat_matrix[begin*ncols +i]);
    }
    return res;
  }
  left = cilk_spawn reduce_max(begin, middle, ncols);
  right = cilk_spawn reduce_max(middle, end, ncols);
  cilk_sync;
  return max(left, right);
}

void fill_histogram(int begin, int end, int ncols) {
  int middle = begin + (end - begin) / 2;
  int i;
  int Self = __cilkrts_get_worker_number();
  if (begin + 1 == end) {
    for (i = 0; i < ncols; i++) {
      histogram[Self][randmat_matrix[begin*ncols +i]]++;
    }
    return;
  }
  cilk_spawn fill_histogram(begin, middle, ncols);
  cilk_spawn fill_histogram(middle, end, ncols);
  cilk_sync;
}

void merge_histogram(int begin, int end) {
  int middle = begin + (end - begin) / 2;
  int i;
  int Cilk_active_size = __cilkrts_get_nworkers();
  if (begin + 1 == end) {
    for (i = 1; i < Cilk_active_size; i++) {
      histogram[0][begin] += histogram[i][begin];
    }
    return;
  }
  cilk_spawn merge_histogram(begin, middle);
  cilk_spawn merge_histogram(middle, end);
  cilk_sync;
}

void fill_mask(int begin, int end, int ncols, int threshold) {
  int middle = begin + (end - begin) / 2;
  int i;
  if (begin + 1 == end) {
    for (i = 0; i < ncols; i++) {
      thresh_mask[begin*ncols +i] = randmat_matrix[begin*ncols +i] >= threshold;
    }
    return;
  }
  cilk_spawn fill_mask(begin, middle, ncols, threshold);
  cilk_spawn fill_mask(middle, end, ncols, threshold);
  cilk_sync;
}

void thresh(int nrows, int ncols, int percent) {
  int i;
  int nmax = 0;
  int count, prefixsum, threshold;

  thresh_mask = (int*) malloc (sizeof(int) * ncols * nrows);

  nmax = cilk_spawn reduce_max(0, nrows, ncols);
  cilk_sync;

  cilk_spawn fill_histogram(0, nrows, ncols);
  cilk_sync;
  cilk_spawn merge_histogram(0, nmax + 1);
  cilk_sync;

  count = (nrows * ncols * percent) / 100;

  prefixsum = 0;
  threshold = nmax;

  for (i = nmax; i >= 0 && prefixsum <= count; i--) {
    prefixsum += histogram[0][i];
    threshold = i;
  }

  cilk_spawn fill_mask(0, nrows, ncols, threshold);
  cilk_sync;
}
