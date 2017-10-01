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

void fill_row(int begin, int ncols, unsigned int seed) {
  int j;
  const int LCG_A = 1664525, LCG_C = 1013904223;
  for (j = 0; j < ncols; j++) {
    seed = LCG_A * seed + LCG_C;
    randmat_matrix[begin*ncols + j] = ((unsigned)seed) % 100;
  }
}

// parallel for on [begin, end), calling f()
void fill_matrix(int begin, int end, int ncols, int seed) {
  int middle = begin + (end - begin) / 2;
  if (begin + 1 == end) {
    cilk_spawn fill_row(begin, ncols, seed + begin);
  } else {
    cilk_spawn fill_matrix(begin, middle, ncols, seed);
    cilk_spawn fill_matrix(middle, end, ncols, seed);
  }
}

void randmat(int nrows, int ncols, int s) {
  // parallel for on rows
  randmat_matrix = (int*) malloc (sizeof(int) * nrows * ncols);
  cilk_spawn fill_matrix(0, nrows, ncols, s);
}
