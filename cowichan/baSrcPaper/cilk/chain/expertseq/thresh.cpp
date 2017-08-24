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
static int histogram[100];

int int_max(int a, int b) {
  return a > b ? a : b;
}

void thresh(int nrows, int ncols, int percent) {
  int i, j;
  int nmax = 0;
  int count, prefixsum, threshold;

  thresh_mask = (int*) malloc (sizeof(int) * ncols * nrows);

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      nmax = int_max(nmax, randmat_matrix[i*ncols + j]);
    }
  }

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      histogram[randmat_matrix[i*ncols + j]]++;
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
      thresh_mask[i*ncols + j] = randmat_matrix[i*ncols + j] >= threshold;
    }
  }
}
