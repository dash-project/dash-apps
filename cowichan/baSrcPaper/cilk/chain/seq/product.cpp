/*
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

void product(int nelts) {
  int i, j;
  double sum;

  product_result = (double*) malloc (sizeof(double) * nelts);

  for (i = 0; i < nelts; i++) {
    sum = 0;
    for (j = 0; j < nelts; j++) {
      sum += outer_matrix[i*nelts + j] * outer_vector[j];
    }
    product_result[i] = sum;
  }
}
