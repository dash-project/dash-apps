/* outer: outer product
 *
 * input:
 *   points: a vector of (x, y) points
 *   nelts: the number of points
 *
 * output:
 *   matrix: a real matrix, whose values are filled with inter-point
 *     distances
 *   vec: a real vector, whose values are filled with origin-to-point
 *     distances
 */
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib> 
#include <cstring>

#include <algorithm>
#include <iostream>
#include <vector>

#include "tbb/tbb.h"

using namespace std;
using namespace tbb;

static int is_bench = 0;
int n_threads = task_scheduler_init::default_num_threads();


static pair<int, int>* points;
static double *matrix;
static double *vec;

typedef blocked_range<size_t> range;

double sqr(double x) {
  return x * x;
}

double distance(const pair<int, int>& x, const pair<int, int>& y) {
  return sqrt(sqr(x.first - y.first) + sqr(x.second - y.second));
}

void outer(int nelts) {
  parallel_for(
    range(0, nelts),
    [&,nelts](range r) {
      for (size_t i = r.begin(); i != r.end(); ++i) {
        double nmax = -1;
        for (int j = 0; j < nelts; j++) {
          if (i != j) {
            matrix[i*nelts + j] = ::distance(points[i], points[j]);
            nmax = max(nmax, matrix[i*nelts + j]);
          }
        }
        matrix[i*nelts + i] = nelts * nmax;
        vec[i] = ::distance(make_pair(0, 0), points[i]);
      }
    });
}

void read_vector_of_points(int nelts) {
  for (int i = 0; i < nelts; i++) {
    cin >> points[i].first >> points[i].second;
  }
}

int main(int argc, char** argv) {
  int nelts;

  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--is_bench")) {
      is_bench = 1;
    } else if (!strcmp(argv[i], "--threads")) {
      sscanf(argv[i + 1], "%d", &n_threads);
      i++;
    }
  }

  task_scheduler_init init(n_threads);

  scanf("%d", &nelts);
  matrix = (double *) malloc (sizeof (double) * nelts * nelts);
  vec = (double *) malloc (sizeof (double) * nelts);
  points = (pair<int, int>*) malloc (sizeof (pair<int, int>) * nelts);

  if (!is_bench) {
    read_vector_of_points(nelts);
  }

  outer(nelts);

  if (!is_bench) {
    printf("%d %d\n", nelts, nelts);
    for (int i = 0; i < nelts; i++) {
      for (int j = 0; j < nelts; j++) {
        printf("%g ", matrix[i*nelts + j]);
      }
      printf("\n");
    }
    printf("\n");

    printf("%d\n", nelts);
    for (int i = 0; i < nelts; i++) {
      printf("%g ", vec[i]);
    }
    printf("\n");
  }

  return 0;
}
