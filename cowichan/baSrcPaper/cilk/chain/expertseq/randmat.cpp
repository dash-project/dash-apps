/*
 * randmat: random number generation
 *
 * input:
 *   nrows, ncols: the number of rows and columns
 *   s: the seed
 *
 * output:
 *   randmat_matrix: an nrows x ncols integer matrix
 */

#include <cilk/cilk.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int *randmat_matrix;

void randmat(int nrows, int ncols, unsigned int s) {
  const int LCG_A = 1664525, LCG_C = 1013904223;
  int i, j;
  randmat_matrix = (int*) malloc (sizeof(int) * nrows * ncols);
  for (i = 0; i < nrows; i++) {
    unsigned seed = s + i;
    for (j = 0; j < ncols; j++) {
      seed = LCG_A * seed + LCG_C;
      randmat_matrix[i*ncols + j] = ((unsigned)seed) % 100;
    }
  }
}
