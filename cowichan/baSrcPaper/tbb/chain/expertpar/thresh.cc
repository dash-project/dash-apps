/* thresh: histogram thresholding
 *
 * input:
 *   randmat_matrix: the integer matrix to be thresholded
 *   nrows, ncols: number of rows and columns
 *   percent: the percentage of cells to retain
 *
 * output:
 *   thresh_mask: a boolean matrix filled with true for cells that are kept
 */
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>

#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"
#include "tbb/parallel_scan.h"
#include "tbb/task_scheduler_init.h"

using namespace std;
using namespace tbb;

extern int is_bench;

extern int *randmat_matrix;
int *thresh_mask;
static int *histogram[100];

typedef tbb::blocked_range<size_t> range;

void thresh(int nrows, int ncols, int percent) {
  int nmax = 0;

  thresh_mask = (int*) malloc (sizeof(int) * nrows * ncols);
  for (int i = 0; i < 100; i++)
    histogram[i] = (int*) malloc (sizeof(int) * nrows);

  nmax = tbb::parallel_reduce(
      range(0, nrows), 0,
      [=](range r, int result)->int {
        for (size_t i = r.begin(); i != r.end(); i++) {
          for (int j = 0; j < ncols; j++) {
            if (is_bench) {
              randmat_matrix[i*ncols + j] = (i * j) % 100;
            }
            result = max(result , (int)randmat_matrix[i*ncols + j]);
          }
        }
        return result;
      },
      [](int x, int y)->int {
        return max(x, y);
      });

  tbb::parallel_for(
      range(0, nrows),
      [=](range r) {
        for (size_t i = r.begin(); i != r.end(); i++) {
          for (int j = 0; j < ncols; j++) {
            histogram[randmat_matrix[i*ncols + j]][i]++;
          }
        }
      });

  tbb::parallel_for(
      range(0, nmax + 1),
      [=](range r) {
        for (size_t j = r.begin(); j != r.end(); j++) {
          for (int i = 1; i < nrows; i++) {
            histogram[j][0] += histogram[j][i];
          }
        }
      });

  int count = (nrows * ncols * percent) / 100;

  int prefixsum = 0;
  int threshold = nmax;

  for (int i = nmax; i >= 0 && prefixsum <= count; i--) {
    prefixsum += histogram[i][0];
    threshold = i;
  }

  tbb::parallel_for(
      range(0, nrows),
      [=](range r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
          for (int j = 0; j < ncols; j++) {
            thresh_mask[i*ncols + j] = randmat_matrix[i*ncols + j] >= threshold;
          }
        }
      });
}
