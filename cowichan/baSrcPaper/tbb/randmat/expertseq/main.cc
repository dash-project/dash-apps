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

static unsigned char* matrix;

int is_bench = 0;

void randmat(int nrows, int ncols, unsigned int s) {
  const int LCG_A = 1664525, LCG_C = 1013904223;
  for (int i = 0; i < nrows; i++) {
    unsigned int seed = s + i;
    for (int j = 0; j < ncols; j++) {
      seed = LCG_A * seed + LCG_C;
      matrix[i*ncols + j] = seed % 100;
    }
  }
}

int main(int argc, char** argv) {
  int nrows, ncols, s;

  if (argc == 2) {
    if (!strcmp(argv[argc - 1], "--is_bench")) {
      is_bench = 1;
    }
  }

  scanf("%d%d%d", &nrows, &ncols, &s);
  matrix = (unsigned char*) malloc (sizeof (unsigned char) * nrows * ncols);
  randmat(nrows, ncols, s);

  if (!is_bench) {
    printf("%d %d\n", nrows, ncols);
    for (int i = 0; i < nrows; i++) {
      for (int j = 0; j < ncols; j++) {
        printf("%d ", matrix[i*ncols + j]);
      }
      printf("\n");
    }
    printf("\n");
  }

  return 0;
}
