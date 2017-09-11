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
#include <cilk/reducer_opadd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <algorithm>

using namespace std;

static int is_bench = 0;

static int *matrix;
static int *mask;
static int *count_per_line;

static pair <int, int> *points;
static pair <int, pair <int, int> > *values;

int reduce_sum (int nrows, int ncols) {
  cilk::reducer_opadd<int> total_sum(0);
  
  cilk_for (int q  = 0; q < nrows; ++q) {
    int tmp_sum = 0;
    for (int i = 0; i < ncols; ++i) {
      // if (is_bench) {
        // mask[q*ncols + i] = ((nrows * i) % (ncols + 1)) == 1;
      // }
      tmp_sum += mask[q * ncols + i];
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
      if (mask[begin*ncols +j] == 1) {
        values[count].first = matrix[begin*ncols +j];
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

  n = reduce_sum(nrows, ncols);

  prefix_sum(nrows + 1);
  fill_values(0, nrows, ncols);

  sort(values, values + n);
  
  // printf("winEl:%i\n",n);

  chunk = n / nelts;
  free (count_per_line);
  points = (pair<int, int>*) malloc (sizeof(pair<int, int>) * nelts);
  cilk_for (int i = 0; i < nelts; i++) {
    index = i * chunk;
    points[i] = values[index].second;
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

inline void FillOnTheFly( int nrows, int ncols, int thresh ){
  int threshInverse = 100 / thresh;
  
  cilk_for( int ix = 0; ix < nrows*ncols; ++ix ){
    if( (ix % threshInverse) == 0 )
    {
      mask[ix] = 1;
    }else{
      mask[ix] = 0;
    }
    
    matrix[ix] = ix % 100;
  }
}

int main(int argc, char *argv[]) {
  int nrows, ncols, nelts, i, thresh, nproc = 0;
  struct timespec start, stop;
  double accum;

  if (argc >= 2) {
    for (int a = 0; a < argc; a++) {
      if (!strcmp(argv[a], "--is_bench")) {
        is_bench = 1;
      } else if (!strcmp(argv[a], "--nproc")) {
        sscanf(argv[++a], "%d", &nproc);
      }
    }
  }

  scanf("%d%d", &nrows, &ncols);
  matrix = (int*) malloc (sizeof(int) * nrows * ncols);
  mask = (int*) malloc (sizeof(int) * nrows * ncols);
  values= (pair<int, pair<int, int> >*) 
  malloc (sizeof(pair<int, pair<int, int> >) * nrows * ncols);

  if (!is_bench) {
    read_matrix(nrows, ncols);
    read_mask(nrows, ncols);
  }

  scanf("%d", &nelts);

  count_per_line = (int*) calloc (nrows+1, sizeof(int));
  
  if(is_bench)
  {
    scanf("%d", &thresh);
    FillOnTheFly( nrows, ncols, thresh );
  }

  if( clock_gettime( CLOCK_MONOTONIC_RAW, &start) == -1 ) {
    perror( "clock gettime error 1" );
    exit( EXIT_FAILURE );
  }
  
  cilk_spawn winnow(nrows, ncols, nelts);
  cilk_sync;
  
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
  fprintf( fp, "Cilk  ,Winnow ,%5i,%5i,%3i,%5i,%2i,%.9lf,isBench:%d\n", nrows, ncols, thresh, nelts, nproc, accum, is_bench );
  fclose ( fp );
  

  if (!is_bench) {
    printf("%d\n", nelts);
    for (i = 0; i < nelts; i++) {
      printf("%d %d\n", points[i].first, points[i].second);
    }
    printf("\n");
  }
  free (matrix);
  free (mask);
  free (values);
  free (points);
  return 0;
}
