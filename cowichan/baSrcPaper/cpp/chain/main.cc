/* chain: chain all problems
 *
 * input:
 *   nelts: the number of elements
 *   randmat_seed: random number generator seed
 *   thresh_percent: percentage of cells to retain
 *   winnow_nelts: the number of points to select
 *
 * output:
 *  result: a real vector, whose values are the result of the final product
 */

#include <cstdio>
#include <iostream>
#include <vector>

using namespace std;

typedef vector<int> VectorInt;
typedef vector<VectorInt> MatrixInt;
typedef vector<double> VectorDouble;
typedef vector<VectorDouble> MatrixDouble;
typedef vector<pair<int, int> > VectorPairInt;

void randmat(int, int, int, MatrixInt*);
void thresh(int, int, const MatrixInt&, int, MatrixInt*);
void winnow(int, int, const MatrixInt&, const MatrixInt&, int,
    VectorPairInt*);
void outer(int, const VectorPairInt&, MatrixDouble*, VectorDouble*);
void product(int, const MatrixDouble&, const VectorDouble&, VectorDouble*);

int main() {
  int nelts, randmat_seed, thresh_percent, winnow_nelts;
  scanf("%d%d%d%d", &nelts, &randmat_seed, &thresh_percent, &winnow_nelts);

  MatrixInt randmat_matrix(nelts, VectorInt(nelts));
  randmat(nelts, nelts, randmat_seed, &randmat_matrix);

  MatrixInt thresh_mask(nelts, VectorInt(nelts));
  thresh(nelts, nelts, randmat_matrix, thresh_percent, &thresh_mask);

  VectorPairInt winnow_points(winnow_nelts);
  winnow(nelts, nelts, randmat_matrix, thresh_mask, winnow_nelts,
      &winnow_points);

  MatrixDouble outer_matrix(winnow_nelts, VectorDouble(winnow_nelts));
  VectorDouble outer_vector(winnow_nelts);
  outer(winnow_nelts, winnow_points, &outer_matrix, &outer_vector);

  VectorDouble product_result(winnow_nelts);
  product(winnow_nelts, outer_matrix, outer_vector, &product_result);

  for (int i = 0; i < winnow_nelts; i++) {
    printf("%g ", product_result[i]);
  }
  printf("\n");

  return 0;
}
