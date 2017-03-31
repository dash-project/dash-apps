/* product: matrix-vector product
 *
 * input:
 *   matrix: a real matrix
 *   vec: a real vector
 *   nelts: the number of elements
 *
 * output:
 *   result: a real vector, whose values are the result of the product
 */
#include <cstdio>
#include <cstring>

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

static int is_bench = 0;

static double *matrix;
static double *vec;
static double *result;

void product(int nelts) {
  for (int i = 0; i < nelts; i++) {
    double sum = 0;
    for (int j = 0; j < nelts; j++) {
      sum = sum + matrix[i*nelts + j] * vec [j];
    }
    result [i] = sum;
  }
}

int main(int argc, char** argv) {
  int nelts;

  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--is_bench")) {
      is_bench = 1;
    }
  }

  scanf("%d", &nelts);

  matrix = (double *) malloc (sizeof(double) * nelts * nelts);
  vec = (double *) malloc (sizeof (double) * nelts);
  result = (double *) malloc (sizeof (double) * nelts);

  if (!is_bench) {
    for (int i = 0; i < nelts; i++) {
      for (int j = 0; j < nelts; j++) {
        scanf("%lf", &matrix[i*nelts + j]);
      }
    }

    for (int i = 0; i < nelts; i++) {
      scanf("%lf", &vec[i]);
    }
  }

  product(nelts);

  if (!is_bench) {
    printf("%d\n", nelts, nelts);
    for (int i = 0; i < nelts; i++) {
      printf("%g ", result[i]);
    }
    printf("\n");
  }

  return 0;
}
