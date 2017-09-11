/* thresh: histogram thresholding
 *
 * input:
 *   matrix: the integer matrix to be thresholded
 *   nrows, ncols: number of rows and columns
 *   percent: the percentage of cells to retain
 *
 * output:
 *   mask: a boolean matrix filled with true for cells that are kept
 */
#include <cstdlib>

#include <algorithm>
#include <vector>

using namespace std;

void thresh(int nrows, int ncols, const vector<vector<int> >& matrix,
    int percent, vector<vector<int> >* mask) {
  int nmax = 0;
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      nmax = max(nmax, matrix[i][j]);
    }
  }

  int* histogram = new int[nmax + 1];

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      histogram[matrix[i][j]]++;
    }
  }

  int count = (nrows * ncols * percent) / 100;

  int prefixsum = 0;
  int threshold = nmax;

  for (int i = nmax; i >= 0 && prefixsum <= count; i--) {
    prefixsum += histogram[i];
    threshold = i;
  }

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      (*mask)[i][j] = matrix[i][j] >= threshold;
    }
  }

  delete[] histogram;
}

