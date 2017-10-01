/*
 * outer: outer product
 *
 * input:
 *   points: a vector of (x, y) points
 *   nelts: the number of points
 *
 * output:
 *   matrix: a real matrix, whose values are filled with inter-point
 *     distances
 *   vector: a real vector, whose values are filled with origin-to-point
 *     distances
 */

#include <cilk/cilk.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static int is_bench = 0;

static double *matrix;
static double *vector;

typedef struct sPoint {
  int i, j;
} Point;

static Point *points;

double sqr(double x) {
  return x * x;
}

double distance(Point a, Point b) {
  return sqrt(sqr(a.i - b.i) + sqr(a.j - b.j));
}

double max(double a, double b) {
  return a > b ? a : b;
}

void outer(int nelts) {
  int i, j;
  double nmax;
  for (i = 0; i < nelts; i++) {
    nmax = -1;
    for (j = 0; j < nelts; j++) {
      if (i != j) {
        matrix[i*nelts + j] = distance(points[i], points[j]);
        nmax = max(nmax, matrix[i*nelts + j]);
      }
      matrix[i*nelts + i] = nmax * nelts;
      vector[i] = distance((Point){0, 0}, points[i]);
    }
  }
}

void read_vector_of_points(int nelts) {
  int i, a, b;
  for (i =  0; i < nelts; i++) {
    scanf("%d %d", &a, &b);
    points[i].i = a;
    points[i].j = b;
  }
}

int main(int argc, char *argv[]) {
  int nelts, i, j;

  if (argc == 2) {
    if (!strcmp(argv[argc - 1], "--is_bench")) {
      is_bench = 1;
    }
  }

  scanf("%d", &nelts);
  matrix = (double*) malloc (sizeof(double) * nelts * nelts);
  vector = (double*) malloc (sizeof(double) * nelts);
  points = (Point*) malloc (sizeof(Point) * nelts);

  if (!is_bench) {
    read_vector_of_points(nelts);
  }

  outer(nelts);

  if (!is_bench) {
    printf("%d %d\n", nelts, nelts);
    for (i = 0; i < nelts; i++) {
      for (j = 0; j < nelts; j++) {
        printf("%g ", matrix[i*nelts + j]);
      }
      printf("\n");
    }
    printf("\n");

    printf("%d\n", nelts);
    for (i = 0; i < nelts; i++) {
      printf("%g ", vector[i]);
    }
    printf("\n");
  }

  return 0;
}
