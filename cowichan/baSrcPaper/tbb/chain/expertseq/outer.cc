/* outer: outer product
 *
 * input:
 *   winnow_points: a vector of (x, y) points
 *   nelts: the number of points
 *
 * output:
 *   outer_matrix: a real matrix, whose values are filled with inter-point
 *     distances
 *   outer_vector: a real vector, whose values are filled with origin-to-point
 *     distances
 */
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib> 
#include <cstring>

#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;

extern pair<int, int> *winnow_points;
double *outer_matrix;
double *outer_vector;

double sqr(double x) {
  return x * x;
}

double distance(const pair<int, int>& x, const pair<int, int>& y) {
  return sqrt(sqr(x.first - y.first) + sqr(x.second - y.second));
}

void outer(int nelts) {
  outer_matrix = (double*) malloc(sizeof(double) * nelts * nelts);
  outer_vector = (double*) malloc(sizeof(double) * nelts);

  for (int i = 0; i < nelts; i++) {
    double nmax = -1;
    for (int j = 0; j < nelts; j++) {
      if (i != j) {
        outer_matrix[i*nelts + j] = ::distance(winnow_points[i], winnow_points[j]);
        nmax = max(nmax, outer_matrix[i*nelts + j]);
      }
    }
    outer_matrix[i*nelts + i] = nelts * nmax;
    outer_vector[i] = ::distance(make_pair(0, 0), winnow_points[i]);
  }
}
