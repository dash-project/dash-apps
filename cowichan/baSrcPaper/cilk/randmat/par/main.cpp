/*
 * randmat: random number generation
 *
 * input:
 *   nrows, ncols: the number of rows and columns
 *   s: the seed
 *
 * output:
 *   matrix: an nrows x ncols integer matrix
 */

// #include <cilk-lib.cilkh>
#include <cilk/cilk.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static int* matrix;

int is_bench = 0;

void fill_row(int begin, int ncols, unsigned int seed) {
  int j;
  const int LCG_A = 1664525, LCG_C = 1013904223;
  for (j = 0; j < ncols; j++) {
    seed = LCG_A * seed + LCG_C;
    matrix[begin*ncols + j] = ((unsigned)seed) % 100;
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
  cilk_spawn fill_matrix(0, nrows, ncols, s);
}

int main(int argc, char *argv[]) {
  int nrows, ncols, s, i, j;

  if (argc == 2) {
    if (!strcmp(argv[argc -1], "--is_bench")) {
      is_bench = 1;
    }
  }

  scanf("%d%d%d", &nrows, &ncols, &s);
  matrix = (int*) malloc (sizeof(int) * nrows * ncols);
  cilk_spawn randmat(nrows, ncols, s);
  cilk_sync;

  if (!is_bench) {
    for (i = 0; i < nrows; i++) {
      for (j = 0; j < ncols; j++) {
        printf("%d ", matrix[i*ncols + j]);
      }
      printf("\n");
    }
    printf("\n");
  }

  return 0;
}
