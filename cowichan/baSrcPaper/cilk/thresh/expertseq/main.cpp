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
static int histogram[100];

int max(int a, int b) {
  return a > b ? a : b;
}

void thresh(int nrows, int ncols, int percent) {
  int i, j;
  int nmax = 0;
  int count, prefixsum, threshold;

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      nmax = max(nmax, matrix[i*ncols + j]);
    }
  }

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      histogram[matrix[i*ncols + j]]++;
    }
  }

  count = (nrows * ncols * percent) / 100;

  prefixsum = 0;
  threshold = nmax;

  for (i = nmax; i >= 0 && prefixsum <= count; i--) {
    prefixsum += histogram[i];
    threshold = i;
  }

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      mask[i*ncols + j] = matrix[i*ncols + j] >= threshold;
    }
  }
}

int main(int argc, char *argv[]) {
  int nrows, ncols, percent, i, j;

  if (argc >= 2) {
    int a;
    for (a = 0; a < argc; a++){
      if (!strcmp(argv[a], "--is_bench")) {
        is_bench = 1;
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

  thresh(nrows, ncols, percent);

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
