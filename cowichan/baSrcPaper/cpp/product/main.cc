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

#include <iostream>
#include <vector>

using namespace std;

void product(int nelts, const vector<vector<double> >& matrix,
    const vector<double>& vec, vector<double>* result) {

  for (int i = 0; i < nelts; i++) {
    double sum = 0;
    for (int j = 0; j < nelts; j++) {
      sum += matrix[i][j] * vec[j];
    }
    (*result)[i] = sum;
  }
}

vector<vector<double> >* matrix;

int main(int argc, char** argv) {
  int nelts;
  scanf("%d", &nelts);

  matrix = new vector<vector<double> >(nelts, vector<double>(nelts));
  vector<double> vec(nelts), result(nelts);

  for (int i = 0; i < nelts; i++) {
    for (int j = 0; j < nelts; j++) {
      cin >> (*matrix)[i][j];
    }
  }

  for (int i = 0; i < nelts; i++) {
    cin >> vec[i];
  }

  product(nelts, *matrix, vec, &result);

  printf("%d\n", nelts, nelts);
  for (int i = 0; i < nelts; i++) {
    printf("%g ", result[i]);
  }
  printf("\n");

  return 0;
}
