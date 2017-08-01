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

#include <cilk/cilk.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static int *matrix;

int is_bench = 0;

void randmat(int nrows, int ncols, int seed) {
  const int LCG_A = 1664525, LCG_C = 1013904223;
  cilk_for (int i = 0; i < nrows; i++) {
    int begin = i;
    int s = seed + i;
    for (int j = 0; j < ncols; j++) {
      s = LCG_A * s + LCG_C;
      matrix[begin*ncols + j] = ((unsigned)s) % 100;
    }
  }
}

int main(int argc, char *argv[]) {
  int nrows, ncols, s, i, j;

  if (argc >= 2) {
    for (int a = 0; a < argc; a++) {
      if (!strcmp(argv[a], "--is_bench")) {
        is_bench = 1;
      }
    }
  }

  scanf("%d%d%d", &nrows, &ncols, &s);
  matrix = (int*) malloc (sizeof(int) * nrows * ncols);
  randmat(nrows, ncols, s);

  //if (!is_bench) {
    for (i = 0; i < nrows; i++) {
      for (j = 0; j < ncols; j++) {
        printf("%d ", matrix[i*ncols + j]);
      }
      printf("\n");
    }
    printf("\n");
  //}

  return 0;
}
