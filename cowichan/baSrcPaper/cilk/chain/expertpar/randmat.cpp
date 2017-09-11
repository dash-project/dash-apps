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

void randmat(int nrows, int ncols, int seed) {
  const int LCG_A = 1664525, LCG_C = 1013904223;
  randmat_matrix = (int*) malloc (sizeof(int) * nrows * ncols);
  cilk_for (int i = 0; i < nrows; i++) {
    int begin = i;
    int s = seed + i;
    for (int j = 0; j < ncols; j++) {
      s = LCG_A * s + LCG_C;
      randmat_matrix[begin*ncols + j] = ((unsigned)s) % 100;
    }
  }
}
