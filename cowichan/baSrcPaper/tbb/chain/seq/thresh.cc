/* thresh: histogram thresholding
 *
 * input:
 *   randmat_matrix: the integer matrix to be thresholded
 *   nrows, ncols: number of rows and columns
 *   percent: the percentage of cells to retain
 *
 * output:
 *   thresh_mask: a boolean matrix filled with true for cells that are kept
 */
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>

using namespace std;

extern int is_bench;
extern int* randmat_matrix;
int* thresh_mask;
static int histogram[100];

void thresh(int nrows, int ncols, int percent) {
  thresh_mask = (int*) malloc (sizeof (int) * ncols * nrows);
  int nmax = 0;
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      if (is_bench) {
        randmat_matrix[i*ncols + j] = (i * j) % 100;
      }
      nmax = max(nmax, randmat_matrix[i*ncols + j]);
    }
  }

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      histogram[randmat_matrix[i*ncols + j]]++;
    }
  }

  int count = (nrows * ncols * percent) / 100;

  int prefixsum = 0;
  int threshold = nmax;

  for (int i = nmax; i >= 0 && prefixsum <= count; i--) {
    prefixsum += histogram[i];
    threshold = i;
  }

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      thresh_mask[i*ncols + j] = randmat_matrix[i*ncols + j] >= threshold;
    }
  }
}
