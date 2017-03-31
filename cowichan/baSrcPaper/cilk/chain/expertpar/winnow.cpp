/*
 * winnow: weighted point selection
 *
 * input:
 *   randmat_matrix: an integer matrix, whose values are used as masses
 *   thresh_mask: a boolean matrix showing which points are eligible for
 *     consideration
 *   nrows, ncols: the number of rows and columns
 *   nelts: the number of points to select
 *
 * output:
 *   winnow_points: a vector of (x, y) points
 */

#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include "common.h"

using namespace std;

extern int is_bench;
extern int *randmat_matrix;
extern int *thresh_mask;
static int *count_per_line;

typedef pair <int, int> point2;

point2 *winnow_points;

static point *values;

int reduce_sum (int nrows, int ncols) {
  cilk::reducer_opadd<int> total_sum(0);
  
  cilk_for (int q  = 0; q < nrows; ++q) {
    int tmp_sum = 0;
    for (int i = 0; i < ncols; ++i) {
      if (is_bench) {
        thresh_mask[q*ncols + i] = ((nrows * i) % (ncols + 1)) == 1;
      }
      tmp_sum += thresh_mask[q * ncols + i];
    } 
    count_per_line[q+1] = tmp_sum;
    total_sum += tmp_sum;
  }
  return total_sum.get_value(); 
}

void scan_update_elements(int begin, int end, int* array, int size) {
  int middle, count;
  if (end - begin <= size) {
    return;
  } else if (end - begin <= 2 * size) {
    array[begin + size] = array[begin] + array[begin + size];
  } else {
    count = (end - begin) / size;
    count /= 2;
    count += count % 2;  // to ensure it is even
    middle = begin + count * size;
    cilk_spawn scan_update_elements(begin, middle, array, size);
    scan_update_elements(middle, end, array, size);
  }
}

// Ladner-Fischer
// parallel scan on [begin, end)
void scan_impl(int begin, int end, int* array, int size) {
  if (end - begin > size) {
    scan_update_elements(begin, end, array, size);
    scan_impl(begin + size, end, array, 2 * size);
    scan_update_elements(begin + size, end, array, size);
  }
}

void scan(int n, int* array) {
  scan_impl(0, n, array, 1);
}

void prefix_sum(int n) {
  scan(n, count_per_line);
}

void fill_values(int begin, int end, int ncols) {
  int middle = begin + (end - begin) / 2;
  int count, j;
  if (begin + 1 == end) {
    count = count_per_line[begin];
    for (j = 0; j < ncols; j++) {
      if (thresh_mask[begin*ncols +j] == 1) {
        values[count].first = randmat_matrix[begin*ncols +j];
        values[count].second.first = begin;
        values[count].second.second = j;
        count++;
      }
    }
    return;
  }
  cilk_spawn fill_values(begin, middle, ncols);
  fill_values(middle, end, ncols);
}

void winnow(int nrows, int ncols, int nelts) {
  int i, n, chunk, index;

  values= (point*) malloc (sizeof(point) * nrows * ncols);
  winnow_points = (point2*) malloc (sizeof(point2) * nelts);
  count_per_line = (int*) calloc (nrows + 1, sizeof(int));

  n = reduce_sum(nrows, ncols);

  prefix_sum(nrows + 1);
  fill_values(0, nrows, ncols);

  sort(values, values + n);

  chunk = n / nelts;

  cilk_for (int i = 0; i < nelts; i++) {
    index = i * chunk;
    winnow_points[i] = values[index].second;
  }

  free (count_per_line);
  free (values);
}
