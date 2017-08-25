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

