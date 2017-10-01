/* thresh: histogram thresholding
 *
 * input:
 *   matrix: the integer matrix to be thresholded
 *   nrows, ncols: number of rows and columns
 *   percent: the percentage of cells to retain
 *
 * output:
 *   mask: a boolean matrix filled with true for cells that are kept
 */
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>

using namespace std;

int is_bench = 0;

static unsigned char* matrix;
static unsigned char* mask;
static int histogram[100];

void thresh(int nrows, int ncols, int percent) {
  int nmax = 0;
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      if (is_bench) {
        matrix[i*ncols + j] = (i * j) % 100;
      }
      nmax = max(nmax, (int)matrix[i*ncols + j]);
    }
  }

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      histogram[matrix[i*ncols + j]]++;
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
      mask[i*ncols + j] = matrix[i*ncols + j] >= threshold;
    }
  }
}

int main(int argc, char** argv) {
  int nrows, ncols, percent;

  if (argc == 2) {
    if (!strcmp(argv[argc - 1], "--is_bench")) {
      is_bench = 1;
    }
  }

  scanf("%d%d", &nrows, &ncols);

  matrix = (unsigned char *) malloc (sizeof (unsigned char) * nrows * ncols);
  mask = (unsigned char *) malloc (sizeof (unsigned char) * nrows * ncols);

  if (!is_bench) {
    for (int i = 0; i < nrows; i++) {
      for (int j = 0; j < ncols; j++) {
        scanf("%hhu", &matrix[i*ncols + j]);
      }
    }
  }

  scanf("%d", &percent);

  thresh(nrows, ncols, percent);

  if (!is_bench) {
    printf("%d %d\n", nrows, ncols);
    for (int i = 0; i < nrows; i++) {
      for (int j = 0; j < ncols; j++) {
        printf("%hhu ", mask[i*ncols + j]);
      }
      printf("\n");
    }
    printf("\n");
  }

  return 0;
}
