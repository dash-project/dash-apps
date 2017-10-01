/*
 * product: matrix-vector product
 *
 * input:
 *   nelts: the number of elements
 *   matrix: a real matrix
 *   vector: a real vector
 *
 * output:
 *   result: a real vector, whose values are the result of the product
 */

#include <cilk/cilk.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int is_bench = 0;

static double *matrix;
static double *vector;
static double *result;


// parallel for on [begin, end)
void fill_result(int begin, int end, int ncols) {
  int middle = begin + (end - begin) / 2;
  double sum = 0;
  int j;
  if (begin + 1 == end) {
    for (j = 0; j < ncols; j++) {
      sum += matrix[begin*ncols + j] * vector[j];
    }
    result[begin] = sum;
    return;
  }
  cilk_spawn fill_result(begin, middle, ncols);
  cilk_spawn fill_result(middle, end, ncols);
}

void product(int nelts) {
  cilk_spawn fill_result(0, nelts, nelts);
}

void read_matrix(int nelts) {
  int i, j;
  for (i = 0; i < nelts; i++) {
    for (j = 0; j < nelts; j++) {
      scanf("%lf", &matrix[i*nelts + j]);
    }
  }
}

void read_vector(int nelts) {
  int i;
  for (i = 0; i < nelts; i++) {
    scanf("%lf", &vector[i]);
  }
}

int main(int argc, char *argv[]) {
  int nelts, i;

  if (argc == 2) {
    if (!strcmp(argv[argc - 1], "--is_bench")) {
      is_bench = 1;
    }
  }

  scanf("%d", &nelts);
  matrix = (double*) malloc (sizeof(double) * nelts * nelts);
  vector = (double*) malloc (sizeof(double) * nelts);
  result = (double*) malloc (sizeof(double) * nelts);

  if (!is_bench) {
    read_matrix(nelts);
    read_vector(nelts);
  }

  cilk_spawn product(nelts);
  cilk_sync;

  if (!is_bench) {
    printf("%d\n", nelts);
    for (i = 0; i < nelts; i++) {
      printf("%g ", result[i]);
    }
    printf("\n");
  }

  return 0;
}
