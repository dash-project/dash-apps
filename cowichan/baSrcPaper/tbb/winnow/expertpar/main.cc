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
#include <cstring>

#include <algorithm>
#include <iostream>
#include <vector>

#include "tbb/tbb.h"

using namespace std;
using namespace tbb;

int is_bench = 0;
int n_threads = task_scheduler_init::default_num_threads();


static int* matrix;
static int* mask;
static pair<int, int> *points;
static pair<int, pair<int, int> > *values;
static int* count_per_line;
static int* total_count;

typedef blocked_range<size_t> range;

class ScanSum {
  int sum;

public:
  ScanSum(): sum(0) {}
  template<typename Tag>
  void operator()(range r, Tag) {
    int res = sum;
    for (size_t i = r.begin(); i != r.end(); ++i) {
      res += count_per_line[i];
      if (Tag::is_final_scan()) {
        total_count[i] = res;
      }
    }
    sum = res;
  }
  ScanSum(ScanSum& other, tbb::split) : sum(0) {}
  void reverse_join(ScanSum& other) { sum += other.sum; }
  void assign(ScanSum& other) { sum = other.sum; }
};

void winnow(int nrows, int ncols, int nelts) {
  int count = 0;

  count = parallel_reduce(
    range(0, nrows), 0,
    [=](range r, int result)->int {
      for (size_t i = r.begin(); i != r.end(); i++) {
        int cur = 0;
        for (int j = 0; j < ncols; j++) {
          if (is_bench) {
            mask[i*ncols + j] = ((i * j) % (ncols + 1)) == 1;
          }
          cur += mask[i*ncols + j];
        }
        result += count_per_line[i + 1] = cur;
      }
      return result;
    },
    [](int x, int y)->int {
      return x + y;
    });

  ScanSum scan_sum;
  tbb::parallel_scan(
      range(0, nrows + 1),
      scan_sum);

  tbb::parallel_for(
      range(0, nrows),
      [=](range r) {
        for (size_t i = r.begin(); i != r.end(); i++) {
          int count = total_count[i];
          for (int j = 0; j < ncols; j++) {
            if (mask[i*ncols + j]) {
              values[count] = (make_pair(matrix[i*ncols + j],
                  make_pair(i, j)));
              count++;
            }
          }
        }
      });


  tbb::parallel_sort (values, values + count);

  size_t n = count;
  size_t chunk = n / nelts;

  for (int i = 0; i < nelts; i++) {
    int index = i * chunk;
    points[i] = values[index].second;
  }
}

void read_matrix(int nrows, int ncols) {
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      int v;
      cin >> v;
      matrix[i*ncols + j] = v;
    }
  }
}

void read_mask(int nrows, int ncols) {
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      int v;
      cin >> v;
      mask[i*ncols + j] = v;
    }
  }
}

int main(int argc, char** argv) {
  int nrows, ncols, nelts;

  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--is_bench")) {
      is_bench = 1;
    } else if (!strcmp(argv[i], "--threads")) {
      sscanf(argv[i + 1], "%d", &n_threads);
      i++;
    }
  }

  task_scheduler_init init(n_threads);

  scanf("%d%d", &nrows, &ncols);
  matrix = (int *) malloc (sizeof (int) * nrows * ncols);
  mask = (int *) malloc (sizeof (int) * nrows * ncols);
  
  total_count = (int *) malloc (sizeof (int) * (nrows + 1));
  memset (total_count, 0, sizeof (int) * (nrows + 1));

  count_per_line = (int *) malloc (sizeof (int) * (nrows + 1));
  memset (count_per_line, 0, sizeof (int) * (nrows + 1));

  if (!is_bench) {
    read_matrix(nrows, ncols);
    read_mask(nrows, ncols);
  }

  scanf("%d", &nelts);
  points = (pair <int, int> *) malloc (sizeof (pair <int, int>) * nelts);
  values = (pair <int, pair <int, int> > *) malloc (sizeof (pair <int, pair <int, int> >) * nrows * ncols);

  winnow(nrows, ncols, nelts);

  if (!is_bench) {
    printf("%d\n", nelts);

    for (int i = 0; i < nelts; i++) {
      printf("%d %d\n", points[i].first, points[i].second);
    }
    printf("\n");
  }

  return 0;
}
