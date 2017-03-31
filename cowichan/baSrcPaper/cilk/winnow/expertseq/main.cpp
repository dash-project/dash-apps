/*
 * winnow: weighted point selection
 *
 * input:
 *   matrix: an integer matrix, whose values are used as masses
 *   mask: a boolean matrix showing which points are eligible for
 *     consideration
 *   nrows, ncols: the number of rows and columns
 *   nelts: the number of points to select
 *
 * output:
 *   points: a vector of (x, y) points
 */

#include <cilk/cilk.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static int is_bench = 0;

static int *matrix;
static int *mask;

typedef struct sPoint {
  int value, i, j;
} Point;

static Point *points;
static Point *values;

int compare(const void* vl, const void* vr) {
  const Point* l = (Point*) vl, *r = (Point *) vr;
  return (l->value - r->value);
}

void winnow(int nrows, int ncols, int nelts) {
  int i, j, n =  0, count = 0, chunk, index;
  for (i = 0; i < nrows; i++) {
    if (is_bench) {
      for (j = 0; j < ncols; j++) {
        mask[i*ncols +j] = ((i * j) % (ncols + 1)) == 1;
      }
    }
    for (j = 0; j < ncols; j++) {
      n += mask[i*ncols +j];
    }
  }

  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      if (mask[i*ncols +j] == 1) {
        values[count].value = matrix[i*ncols +j];
        values[count].i = i;
        values[count].j = j;
        count++;
      }
    }
  }
  free (mask);
  free (matrix);

  qsort(values, n, sizeof(*values), compare);

  chunk = n / nelts;

  points = (Point*) malloc (sizeof(Point) * nelts);

  for (i = 0; i < nelts; i++) {
    index = i * chunk;
    points[i] = values[index];
  }
  
}

void read_matrix(int nrows, int ncols) {
  int i, j;
  for (i =  0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      scanf("%hhu", &matrix[i*ncols +j]);
    }
  }
}

void read_mask(int nrows, int ncols) {
  int i, j;
  for (i =  0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
      scanf("%hhu", &mask[i*ncols +j]);
    }
  }
}

int main(int argc, char *argv[]) {
  int nrows, ncols, nelts, i;

  if (argc >= 2) {
    int a;
    for (a = 0; a < argc; a++){
      if (!strcmp(argv[a], "--is_bench")) {
        is_bench = 1;
      }
    }
  }

  scanf("%d%d", &nrows, &ncols);

  matrix = (int*) malloc (sizeof(int) * nrows * ncols);
  mask = (int*) malloc (sizeof(int) * nrows * ncols);
  values= (Point*) malloc (sizeof(Point) * nrows * ncols);

  if (!is_bench) {
    read_matrix(nrows, ncols);
    read_mask(nrows, ncols);
  }

  scanf("%d", &nelts);

  winnow(nrows, ncols, nelts);

  if (!is_bench) {
    printf("%d\n", nelts);
    for (i = 0; i < nelts; i++) {
      printf("%d %d\n", points[i].i, points[i].j);
    }
    printf("\n");
  }

  return 0;
}
