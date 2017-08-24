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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int is_bench = 0;

static int *matrix;
static int *mask;
static int histogram[256][100];

int max(int a, int b) {
  return a > b ? a : b;
}

// parallel max reduce on [begin, end)
int reduce_max(int begin, int end, int ncols) {
  int middle = begin + (end - begin) / 2;
  int left, right, res, i;
  if (begin + 1 == end) {
    res = matrix[begin*ncols];
    for (i = 1; i < ncols; i++) {
      res = max(res, matrix[begin*ncols + i]);
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
      histogram[Self][matrix[begin*ncols + i]]++;
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
      mask[begin*ncols + i] = matrix[begin*ncols + i] >= threshold;
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

int main(int argc, char *argv[]) {
  int nrows, ncols, percent, i, j;

  if (argc == 2) {
    if (!strcmp(argv[argc - 1], "--is_bench")) {
      is_bench = 1;
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

  cilk_spawn thresh(nrows, ncols, percent);
  cilk_sync;

  if (!is_bench) {
    printf("%d %d\n", nrows, ncols);
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
