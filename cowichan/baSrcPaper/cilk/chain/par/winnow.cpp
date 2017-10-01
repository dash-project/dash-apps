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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

extern int is_bench;
extern int *randmat_matrix;
extern int *thresh_mask;
static int *count_per_line;

typedef struct sPoint {
  int value, i, j;
} Point;

Point *winnow_points;
static Point *values;

int compare(const void* vl, const void* vr) {
  const Point* l = (Point *) vl, *r = (Point *) vr;
  return (l->value - r->value);
}

int reduce_sum(int begin, int end, int ncols) {
  int middle = begin + (end - begin) / 2;
  int left, right, res, i;
  if (begin + 1 == end) {
    if (is_bench) {
      for (i = 0; i < ncols; i++) {
        thresh_mask[begin*ncols + i] = ((begin * i) % (ncols + 1)) == 1;
      }
    }
    res = thresh_mask[begin*ncols + 0];
    for (i = 1; i < ncols; i++) {
      res += thresh_mask[begin*ncols + i];
    }
    return count_per_line[begin + 1] = res;
  }
  left = cilk_spawn reduce_sum(begin, middle, ncols);
  right = cilk_spawn reduce_sum(middle, end, ncols);
  cilk_sync;
  return left + right;
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
    cilk_spawn scan_update_elements(middle, end, array, size);
  }
}

// Ladner-Fischer
// parallel scan on [begin, end)
void scan_impl(int begin, int end, int* array, int size) {
  if (end - begin > size) {
    cilk_spawn scan_update_elements(begin, end, array, size);
    cilk_sync;
    cilk_spawn scan_impl(begin + size, end, array, 2 * size);
    cilk_sync;
    cilk_spawn scan_update_elements(begin + size, end, array, size);
    cilk_sync;
  }
}

void scan(int n, int* array) {
  cilk_spawn scan_impl(0, n, array, 1);
}

void prefix_sum(int n) {
  cilk_spawn scan(n, count_per_line);
}

void fill_values(int begin, int end, int ncols) {
  int middle = begin + (end - begin) / 2;
  int count, j;
  if (begin + 1 == end) {
    count = count_per_line[begin];
    for (j = 0; j < ncols; j++) {
      if (thresh_mask[begin*ncols + j] == 1) {
        values[count].value = randmat_matrix[begin*ncols + j];
        values[count].i = begin;
        values[count].j = j;
        count++;
      }
    }
    return;
  }
  cilk_spawn fill_values(begin, middle, ncols);
  cilk_spawn fill_values(middle, end, ncols);
}

void winnow(int nrows, int ncols, int nelts) {
  int i, n =  0, chunk, index;

  values = (Point*) malloc (sizeof(Point) * nrows * ncols);
  winnow_points = (Point*) malloc (sizeof(Point) * nelts);
  count_per_line = (int*) malloc (sizeof(int) * (nrows + 1));

  n = cilk_spawn reduce_sum(0, nrows, ncols);
  cilk_sync;

  cilk_spawn prefix_sum(nrows + 1);
  cilk_sync;

  cilk_spawn fill_values(0, nrows, ncols);
  cilk_sync;

  qsort(values, n, sizeof(*values), compare);

  chunk = n / nelts;

  for (i = 0; i < nelts; i++) {
    index = i * chunk;
    winnow_points[i] = values[index];
  }

  free (count_per_line);
  free (values);
}
