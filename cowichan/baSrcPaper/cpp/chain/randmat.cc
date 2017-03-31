/* randmat: random number generator
 *
 * input:
 *   nrows, ncols: number of rows and columns
 *   s: random number generation seed
 *
 * output:
 *   matrix: random nrows x ncols integer matrix
 */
#include <cstdlib>

#include <vector>

using namespace std;

#define MAX_NUMBER 1000

void randmat(int nrows, int ncols, int s, vector<vector<int> >* matrix) {
  srand(s);
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      (*matrix)[i][j] = rand() % MAX_NUMBER;
    }
  }
}

