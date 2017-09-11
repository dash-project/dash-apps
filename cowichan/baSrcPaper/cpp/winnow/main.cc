/* winnow: weighted point selection
 *
 * input:
 *   matrix: an integer matrix, whose values are used as masses
 *   mask: a boolean matrix, showing which points are eligible for
 *     consideration
 *   nrows, ncols: number of rows and columns
 *   nelts: the number of points to select
 *
 * output:
 *   points: a vector of (x, y) points
 */
#include <cassert>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;

void winnow(int nrows, int ncols, const vector<vector<int> >& matrix,
    const vector<vector<int> >& mask, int nelts,
    vector<pair<int, int> >* points) {
  vector<pair<int, pair<int, int> > > values;
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      if (mask[i][j]) {
        values.push_back(make_pair(matrix[i][j], make_pair(i, j)));
      }
    }
  }
  sort(values.begin(), values.end());

  size_t n = values.size();
  size_t chunk = n / nelts;

  for (int i = 0; i < nelts; i++) {
    int index = i * chunk;
    (*points)[i] = values[index].second;
  }

}

void read_matrix(int nrows, int ncols, vector<vector<int> >* matrix) {
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      cin >> (*matrix)[i][j];
    }
  }
}

vector<vector<int> >* matrix;
vector<vector<int> >* mask;

int main(int argc, char** argv) {
  int nrows, ncols, nelts;

  scanf("%d%d", &nrows, &ncols);

  matrix = new vector<vector <int> >(nrows, vector<int>(ncols));
  mask = new vector<vector<int> >(nrows, vector<int>(ncols));

  read_matrix(nrows, ncols, matrix);
  read_matrix(nrows, ncols, mask);

  scanf("%d", &nelts);

  vector<pair<int, int> > points(nelts);

  winnow(nrows, ncols, *matrix, *mask, nelts, &points);

  printf("%d\n", nelts);

  for (int i = 0; i < nelts; i++) {
    printf("%d %d\n", points[i].first, points[i].second);
  }
  printf("\n");

  return 0;
}
