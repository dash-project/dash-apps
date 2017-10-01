#ifndef __COMMON_H__
#define __COMMON_H__
#include <algorithm>

using namespace std; 

typedef pair <int, pair <int, int> > point;

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
  point* values;
  int* count;
  point* points;
  int chunk;
  int nelts;
  double** matrix_double;
  double* vector;
  double* result;
} SplitState;

typedef int (*GenericReduceFunction)(int index, ReduceState state);
typedef void (*SplitFunction)(int index, SplitState state);

void split(int begin, int end, SplitFunction f, SplitState state);

int generic_reduce(int begin, int end, Operator op,
    GenericReduceFunction f, ReduceState state);

int reduce2d_identity_filter(int col, ReduceState state);

int reduce2d_start_cols(int row, ReduceState state);

int reduce2d_with_filter(int nrows, int ncols, int** matrix,
    Operator op, GenericReduceFunction filter, int extra);

int reduce2d(int nrows, int ncols, int** matrix, Operator op);

int filter_matrix(int col, ReduceState state);

int sum(int a, int b);

// Prototypes for the individual pieces of the pipeline
void randmat(int, int, int);
void thresh(int, int, int);
void winnow(int, int, int);
void outer(int);
void product(int);


#endif
