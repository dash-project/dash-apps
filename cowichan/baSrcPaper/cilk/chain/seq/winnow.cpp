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

typedef struct sPoint {
  int value, i, j;
} Point;

Point *winnow_points;
static Point *values;

int compare(const void* vl, const void* vr) {
  const Point* l = (Point*) vl, *r = (Point *) vr;
  return (l->value - r->value);
}

void winnow(int nrows, int ncols, int nelts) {
  int i, j, n =  0, count = 0, chunk, index;
  values = (Point*) malloc (sizeof(Point) * nrows * ncols);
  winnow_points = (Point*) malloc (sizeof(Point) * nrows * ncols);

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      n += thresh_mask[i*ncols +j];
    }
  }

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      if (thresh_mask[i*ncols +j] == 1) {
        values[count].value = randmat_matrix[i*ncols +j];
        values[count].i = i;
        values[count].j = j;
        count++;
      }
    }
  }

  qsort(values, n, sizeof(*values), compare);

  chunk = n / nelts;

  for (i = 0; i < nelts; i++) {
    index = i * chunk;
    winnow_points[i] = values[index];
  }

  free (values);
}
