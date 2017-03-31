/*
 * chain: chain all problems
 *
 * input:
 *   nelts: the number of elements
 *   randmat_seed: random number generator of cells to retain
 *   thresh_percent: percentage of cells to retain
 *   winnow_nelts: the number of points to select
 *
 * output:
 *   result: a real vector, whose values are the result of the final product
 */
#include <cilk/cilk.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <cilk/cilk.h>

int is_bench = 0;

void randmat(int, int, int);
void thresh(int, int, int);
void winnow(int, int, int);
void outer(int);
void product(int);

extern double *product_result;

int main(int argc, char** argv) {
  int nelts, randmat_seed, thresh_percent, winnow_nelts, i;

  if (argc == 2) {
    if (!strcmp(argv[argc - 1], "--is_bench")) {
      is_bench = 1;
    }
  }

  scanf("%d%d%d%d", &nelts, &randmat_seed, &thresh_percent, &winnow_nelts);

  cilk_spawn randmat(nelts, nelts, randmat_seed); cilk_sync;
  cilk_spawn thresh(nelts, nelts, thresh_percent); cilk_sync;
  cilk_spawn winnow(nelts, nelts, winnow_nelts); cilk_sync;
  cilk_spawn outer(winnow_nelts); cilk_sync;
  cilk_spawn product(winnow_nelts); cilk_sync;

  if (!is_bench) {
    for (i = 0; i < winnow_nelts; i++) {
      printf("%g ", product_result[i]);
    }
    printf("\n");
  }

  return 0;
}

