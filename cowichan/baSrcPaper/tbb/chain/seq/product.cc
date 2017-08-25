/* product: matrix-vector product
 *
 * input:
 *   outer_matrix: a real matrix
 *   outer_vector: a real vector
 *   nelts: the number of elements
 *
 * output:
 *   product_result: a real vector, whose values are the result of the product
 */
#include <cstdio>
#include <cstring>

#include <iostream>
#include <vector>

using namespace std;

extern double* outer_matrix;
extern double* outer_vector;
double *product_result;

void product(int nelts) {
  product_result = (double*) malloc (sizeof(double) * nelts);
  for (int i = 0; i < nelts; i++) {
    double sum = 0;
    for (int j = 0; j < nelts; j++) {
      sum += outer_matrix[i*nelts + j] * outer_vector[j];
    }
    product_result[i] = sum;
  }
}
