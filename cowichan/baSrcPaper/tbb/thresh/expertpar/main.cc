/* thresh: histogram thresholding
 *
 * input:
 *   matrix: the integer matrix to be thresholded
 *   nrows, ncols: number of rows and columns
 *   percent: the percentage of cells to retain
 *
 * output:
 *   mask: a boolean matrix filled with true for cells that are kept
 */
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <time.h>
#include <stdio.h>

#include <algorithm>

#include "tbb/tbb.h"

using namespace std;
using namespace tbb;

int is_bench = 0;
int n_threads = task_scheduler_init::default_num_threads();

static unsigned char *matrix;
static unsigned char *mask;

struct view {
  int h[100];
  view() {fill_n (h, 100, 0);}
};

typedef combinable<view> histogram_type;

typedef tbb::blocked_range<size_t> range;

void thresh(int nrows, int ncols, int percent) {
  int nmax = 0;
  histogram_type histogram;

  nmax = tbb::parallel_reduce(
      range(0, nrows), 0,
      [=,&histogram](range r, int result)->int {
        view& v =  histogram.local ();
        for (size_t i = r.begin(); i != r.end(); i++) {
          for (int j = 0; j < ncols; j++) {
            int val;
            if (is_bench) {
              matrix[i*ncols + j] = (i * j) % 100;
            }
            val = (int)matrix[i*ncols + j];

            result = max(result, val);
            v.h[val]++;
          }
        }
        return result;
      },
      [](int x, int y)->int {
        return max(x, y);
      });

  view v;
  histogram.combine_each( [=, &v] (const view& x) {
      for (int i = 0; i <= nmax; ++i)
        v.h[i] += x.h[i];
    });

  int count = (static_cast<size_t>(nrows) * ncols * percent) / 100;
  // printf("thresh:%i\n",count);

  int prefixsum = 0;
  int threshold = nmax;

  for (int i = nmax; i >= 0 && prefixsum <= count; i--) {
    prefixsum += v.h[i];
    threshold = i;
  }

  tbb::parallel_for(
      range(0, nrows),
      [=](range r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
          for (int j = 0; j < ncols; j++) {
            mask[i*ncols + j] = matrix[i*ncols + j] >= threshold;
          }
        }
      });
}

int main(int argc, char** argv) {
  int nrows, ncols, percent;
  struct timespec start, stop;
  double accum;

  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--is_bench")) {
      is_bench = 1;
    } else if (!strcmp(argv[i], "--threads")) {
      sscanf(argv[i + 1], "%d", &n_threads);
      i++;
    }
  }

  task_scheduler_init init(n_threads);

  scanf("%d%d", &nrows, &ncols);
  
  matrix = (unsigned char*) malloc (sizeof (unsigned char) * nrows * ncols);
  mask = (unsigned char*) malloc (sizeof (unsigned char) * nrows * ncols);

  if (!is_bench) {
    for (int i = 0; i < nrows; i++) {
      for (int j = 0; j < ncols; j++) {
        scanf("%hhu", &matrix[i*ncols + j]);
      }
    }
  }

  scanf("%d", &percent);
  
  if( clock_gettime( CLOCK_MONOTONIC_RAW, &start) == -1 ) {
    perror( "clock gettime error 1" );
    exit( EXIT_FAILURE );
  }

  thresh(nrows, ncols, percent);
  
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
  fprintf( fp, "TBB   ,Thresh ,%5u,%5u,%3u,     ,%2u,%.9lf, isBench:%d\n", nrows, ncols, percent, n_threads, accum, is_bench );
  fclose ( fp );
  

  if (!is_bench) {
    // printf("%d %d\n", nrows, ncols);
    for (int i = 0; i < nrows; i++) {
      for (int j = 0; j < ncols; j++) {
        printf("%hhu ", mask[i*ncols + j]);
      }
      printf("\n");
    }
    printf("\n");
  }

  return 0;
}
