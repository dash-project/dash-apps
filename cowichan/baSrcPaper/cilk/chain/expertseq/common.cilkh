#ifndef __COMMON_H__
#define __COMMON_H__

typedef struct sPoint {
  int value, i, j;
} Point;

typedef int (*Operator)(int, int);
typedef void* ReduceFilter;

typedef struct sReduceState {
  int extra, ncols;
  int row;
  int** matrix;
  Operator op;
  ReduceFilter filter;
} ReduceState;

typedef struct sSplitState {
  int* histogram;
  int nrows, ncols;
  int** matrix_int;
  int row;
  int** mask;
  int threshold;
  Point* values;
  int* count;
  Point* points;
  int chunk;
  int nelts;
  double** matrix_double;
  double* vector;
  double* result;
} SplitState;

typedef cilk int (*GenericReduceFunction)(int index, ReduceState state);
typedef cilk void (*SplitFunction)(int index, SplitState state);

cilk void split(int begin, int end, SplitFunction f, SplitState state);

cilk int generic_reduce(int begin, int end, Operator op,
    GenericReduceFunction f, ReduceState state);

cilk int reduce2d_identity_filter(int col, ReduceState state);

cilk int reduce2d_start_cols(int row, ReduceState state);

cilk int reduce2d_with_filter(int nrows, int ncols, int** matrix,
    Operator op, GenericReduceFunction filter, int extra);

cilk int reduce2d(int nrows, int ncols, int** matrix, Operator op);

cilk int filter_matrix(int col, ReduceState state);

int sum(int a, int b);

#endif
