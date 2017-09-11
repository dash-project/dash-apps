/* randmat: random number generator
 *
 * input:
 *   nrows, ncols: number of rows and columns
 *   s: random number generation seed
 *
 * output:
 *   randmat_matrix: random nrows x ncols integer matrix
 */
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/task_scheduler_init.h"

using namespace tbb;

typedef blocked_range<size_t> range;

int *randmat_matrix;

void randmat(int nrows, int ncols, unsigned int s) {
  const int LCG_A = 1664525, LCG_C = 1013904223;
  randmat_matrix = (int*) malloc (sizeof (int) * ncols * nrows);
  parallel_for(
    range(0, nrows),
    [=](range r) {
      for (size_t i = r.begin(); i != r.end(); ++i) {
        unsigned int seed = s + i;
        for (int j = 0; j < ncols; j++) {
          seed = LCG_A * seed + LCG_C;
          randmat_matrix[i*ncols + j] = seed % 100;
        }
      }
  });
}
