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
#include <cstdlib>

#include <algorithm>
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

