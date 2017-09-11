/*
 * outer: outer product
 *
 * input:
 *   winnow_points: a vector of (x, y) points
 *   nelts: the number of points
 *
 * output:
 *   outer_matrix: a real matrix, whose values are filled with inter-point
 *     distances
 *   outer_vector: a real vector, whose values are filled with origin-to-point
 *     distances
 */

#include <cilk/cilk.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

double *outer_matrix;
double *outer_vector;

typedef struct sPoint {
  int v, i, j;
} Point;

extern Point *winnow_points;

double sqr(double x) {
  return x * x;
}

double distance(Point a, Point b) {
  return sqrt(sqr(a.i - b.i) + sqr(a.j - b.j));
}

double double_max(double a, double b) {
  return a > b ? a : b;
}

void outer(int nelts) {
  int i, j;
  double nmax;

  outer_matrix = (double*) malloc (sizeof(double) * nelts * nelts);
  outer_vector = (double*) malloc (sizeof(double) * nelts);

  for (i = 0; i < nelts; i++) {
    nmax = -1;
    for (j = 0; j < nelts; j++) {
      if (i != j) {
        outer_matrix[i*nelts + j] = distance(winnow_points[i], winnow_points[j]);
        nmax = double_max(nmax, outer_matrix[i*nelts + j]);
      }
      outer_matrix[i*nelts + i] = nmax * nelts;
      outer_vector[i] = distance((Point){0, 0, 0}, winnow_points[i]);
    }
  }
}
