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
#include <time.h>

#include <algorithm>

using namespace std;

static int is_bench = 0;

static double *matrix;
static double *vector;

typedef std::pair<int, int> point;
static point *points;

double sqr(double x) {
  return x * x;
}

double distance(point a, point b) {
  return sqrt(sqr(a.first - b.first) + sqr(a.second - b.second));
}

// parallel for on [begin, end)
void outer (int nelts) {
  cilk_for (int i = 0; i < nelts; ++i) {
    double nmax = 0;
    for (int j = 0; j < nelts; ++j) {
      if (i != j) {
        matrix [i*nelts + j] = ::distance (points [i], points [j]);
        nmax = max (nmax, matrix [i*nelts + j]);
      }      
    }
    matrix [i*nelts + i] = nmax * nelts;
    vector [i] = ::distance (make_pair (0,0), points [i]);
  }
}


void read_vector_of_points(int nelts) {
  int i, a, b;
  for (i =  0; i < nelts; i++) {
    scanf("%d %d", &a, &b);
    points[i] = make_pair (a, b);
  }
}

int main(int argc, char *argv[]) {
  int nelts, i, j, nproc = 0;
  struct timespec start, stop;
  double accum;

  if (argc >= 2) {
    for (int a = 0; a < argc; a++){
      if (!strcmp(argv[a], "--is_bench")) {
        is_bench = 1;
      } else if (!strcmp(argv[a], "--nproc")) {
        sscanf(argv[++a], "%d", &nproc);
      }
    }
  }

  scanf("%d", &nelts);
  matrix = (double*) malloc (sizeof(double) * nelts * nelts);
  vector = (double*) malloc (sizeof(double) * nelts);
  points = (point*) malloc (sizeof(point) * nelts);


  if (!is_bench) {
    read_vector_of_points(nelts);
  }
  
  if( clock_gettime( CLOCK_MONOTONIC_RAW, &start) == -1 ) {
    perror( "clock gettime error 1" );
    exit( EXIT_FAILURE );
  }

  outer(nelts);
  
  if( clock_gettime( CLOCK_MONOTONIC_RAW, &stop) == -1 ) {
    perror( "clock gettime error 2" );
    exit( EXIT_FAILURE );
  }
  
  
  accum = ( stop.tv_sec - start.tv_sec ) + ( stop.tv_nsec - start.tv_nsec ) / 1e9;
  

  FILE* fp = fopen("./measurements.txt", "a");
  
  if( !fp ) {
      perror("File opening for benchmark results failed");
      return EXIT_FAILURE;
  }
  // Lang, Problem, rows, cols, thresh, winnow_nelts, jobs, time
  fprintf( fp, "Cilk  ,Outer  ,     ,     ,   ,%5u,%2u,%.9lf,isBench:%d\n", nelts, nproc, accum, is_bench );
  fclose ( fp );
  

  if (!is_bench) {
    printf("%d\n", nelts);
    for (i = 0; i < nelts; i++) {
      for (j = 0; j < nelts; j++) {
        printf("%.4f ", matrix[i*nelts + j]);
      }
      printf("\n");
    }
    printf("\n");

    for (i = 0; i < nelts; i++) {
      printf("%.4f ", vector[i]);
    }
    printf("\n");
  }

  return 0;
}
