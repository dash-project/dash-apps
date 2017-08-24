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

void randmat(int nrows, int ncols, unsigned int s) {
  int i, j;
  const int LCG_A = 1664525, LCG_C = 1013904223;
  unsigned int seed;
  for (i = 0; i < nrows; i++) {
    seed = s + i;
    for (j = 0; j < ncols; j++) {
      seed = LCG_A * seed + LCG_C;
      matrix[i*ncols + j] = ((unsigned)seed) % 100;
    }
  }
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
  randmat(nrows, ncols, s);

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
