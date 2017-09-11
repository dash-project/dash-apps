/* randmat: random number generator
 *
 * input:
 *   nrows, ncols: number of rows and columns
 *   s: random number generation seed
 *
 * output:
 *   matrix: random nrows x ncols integer matrix
 */
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

int* randmat_matrix;

void randmat(int nrows, int ncols, unsigned int s) {
  randmat_matrix = (int*) malloc (sizeof(int) * ncols * nrows);
  const int LCG_A = 1664525, LCG_C = 1013904223;
  for (int i = 0; i < nrows; i++) {
    unsigned int seed = s + i;
    for (int j = 0; j < ncols; j++) {
      seed = LCG_A * seed + LCG_C;
      randmat_matrix[i*ncols + j] = seed % 100;
    }
  }
}
