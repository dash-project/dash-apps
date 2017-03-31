/* outer: outer product
 *
 * input:
 *   points: a vector of (x, y) points
 *   nelts: the number of points
 *
 * output:
 *   matrix: a real matrix, whose values are filled with inter-point
 *     distances
 *   vec: a real vector, whose values are filled with origin-to-point
 *     distances
 */
#include <cmath>
#include <cstdlib>

#include <algorithm>
#include <vector>

using namespace std;

double sqr(double x) {
  return x * x;
}

double distance(const pair<int, int>& x, const pair<int, int>& y) {
  return sqrt(sqr(x.first - y.first) + sqr(x.second - y.second));
}

void outer(int nelts, const vector<pair<int, int> > & points,
    vector<vector<double> >* matrix, vector<double>* vec) {
  for (int i = 0; i < nelts; i++) {
    double nmax = -1;
    for (int j = 0; j < nelts; j++) {
      if (i != j) {
        (*matrix)[i][j] = ::distance(points[i], points[j]);
        nmax = max(nmax, (*matrix)[i][j]);
      }
    }
    (*matrix)[i][i] = nelts * nmax;
    (*vec)[i] = ::distance(make_pair(0, 0), points[i]);
  }
}

