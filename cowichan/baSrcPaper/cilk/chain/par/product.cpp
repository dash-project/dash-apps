/*
 * product: matrix-vector product
 *
 * input:
 *   nelts: the number of elements
 *   outer_matrix: a real matrix
 *   outer_vector: a real vector
 *
 * output:
 *   product_result: a real vector, whose values are the result of the product
 */

#include <cilk/cilk.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern double *outer_matrix;
extern double *outer_vector;
double *product_result;

// parallel for on [begin, end)
void fill_result(int begin, int end, int ncols) {
  int middle = begin + (end - begin) / 2;
  double sum = 0;
  int j;
  if (begin + 1 == end) {
    for (j = 0; j < ncols; j++) {
      sum += outer_matrix[begin*ncols + j] * outer_vector[j];
    }
    product_result[begin] = sum;
    return;
  }
  cilk_spawn fill_result(begin, middle, ncols);
  cilk_spawn fill_result(middle, end, ncols);
}

void product(int nelts) {
  product_result = (double*) malloc (sizeof (double) * nelts);
  cilk_spawn fill_result(0, nelts, nelts);
}
