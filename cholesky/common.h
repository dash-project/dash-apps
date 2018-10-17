#ifndef CHOLESKY_COMMON_H
#define CHOLESKY_COMMON_H

#include <mkl.h>
#include "MatrixBlock.h"
#include "ExtraeInstrumentation.h"

typedef dash::util::Timer<dash::util::TimeMeasure::Clock> Timer;

extern "C" {
void dgemm_ (const char *transa, const char *transb, int *l, int *n, int *m, double *alpha,
             const void *a, int *lda, void *b, int *ldb, double *beta, void *c, int *ldc);
void dtrsm_ (char *side, char *uplo, char *transa, char *diag, int *m, int *n, double *alpha,
             double *a, int *lda, double *b, int *ldb);
void dsyrk_ (char *uplo, char *trans, int *n, int *k, double *alpha, double *a, int *lda,
             double *beta, double *c, int *ldc);
//void dpotrf_(const char* uplo, int& n, float* a, int* lda, int* info);
}

static
void potrf(double * const A, int ts, int ld)
{
   static int INFO;
   static char L = 'L';
   EXTRAE_ENTER(EVENT_POTRF);
   dpotrf_(&L, &ts, A, &ld, &INFO);
   EXTRAE_EXIT(EVENT_POTRF);
}

static
void trsm(double *A, double *B, int ts, int ld)
{
   static char LO = 'L', TR = 'T', NU = 'N', RI = 'R';
   static double DONE = 1.0;
   EXTRAE_ENTER(EVENT_TRSM);
   dtrsm_(&RI, &LO, &TR, &NU, &ts, &ts, &DONE, A, &ld, B, &ld );
   EXTRAE_EXIT(EVENT_TRSM);
}

static
void syrk(double *A, double *B, int ts, int ld)
{
   static char LO = 'L', NT = 'N';
   static double DONE = 1.0, DMONE = -1.0;
   EXTRAE_ENTER(EVENT_SYRK);
   dsyrk_(&LO, &NT, &ts, &ts, &DMONE, A, &ld, &DONE, B, &ld );
   EXTRAE_EXIT(EVENT_SYRK);
}

static
void gemm(double *A, double *B, double *C, int ts, int ld)
{
   static const char TR = 'T', NT = 'N';
   static double DONE = 1.0, DMONE = -1.0;
   EXTRAE_ENTER(EVENT_GEMM);
   dgemm_(&NT, &TR, &ts, &ts, &ts, &DMONE, A, &ld, B, &ld, &DONE, C, &ld);
   EXTRAE_EXIT(EVENT_GEMM);
}

template<typename MatrixT>
static
void
compute_single(MatrixT& matrix, size_t block_size){

  using LocalBlockCache = MatrixBlock<MatrixT>;

  const size_t num_blocks = matrix.pattern().blockspec().extent(0);
  /**
   * Algorithm taken from
   * https://pm.bsc.es/ompss-docs/examples/01-examples/cholesky
   */

  if (dash::myid() == 0){
    // iterate over column of blocks
    for (size_t k = 0; k < num_blocks; ++k) {

      if (dash::myid() == 0)
        std::cout << "Processing column " << k << std::endl;

      LocalBlockCache block_k(matrix, k, k);

      /**
      * Compute cholesky factorization of block on diagonal
      */
      potrf(block_k.lbegin(), block_size, block_size);
      block_k.store();

      /**
      * Solve the triangular equation system in the block
      */
      for (int i = k+1; i < num_blocks; ++i) {
        LocalBlockCache block_b(matrix, k, i);
        trsm(block_k.lbegin(),
             block_b.lbegin(), block_size, block_size);
        block_b.store();
      }

      // walk to the right
      for (size_t i = k+1; i < num_blocks; ++i) {
        // run down to the diagonal
        LocalBlockCache block_a(matrix, k, i);
        for (size_t j = k+1; j < i; ++j) {
          LocalBlockCache block_c(matrix, j, i);
          LocalBlockCache block_b(matrix, k, j);
          // A[k,i] = A[k,i] - A[k,j] * (A[j,i])^t
          gemm(block_a.lbegin(),
               block_b.lbegin(),
               block_c.lbegin(), block_size, block_size);
          block_c.store();
        }
      }

      // update diagonal blocks
      for (size_t i = k+1; i < num_blocks; ++i) {
        LocalBlockCache block_i(matrix, i, i);
        LocalBlockCache block_a(matrix, k, i);
        // A[j,j] = A[j,j] - A[j,i] * (A[j,i])^t
        syrk(block_a.lbegin(),
             block_i.lbegin(), block_size, block_size);
        block_i.store();
      }
    }
  }

  dash::barrier();
}


template<typename MatrixT>
static
void print_matrix(MatrixT &matrix)
{
  if (matrix.extent(0) > 100) return;

  using value_t = typename MatrixT::value_type;

  std::cout << std::setw(6);

  if (dash::myid() == 0) {
    for (size_t i = 0; i < matrix.extent(0); ++i) {
      for (size_t j = 0; j < matrix.extent(1); ++j) {
        std::cout << std::setw(5) << std::setprecision(3)
                  << static_cast<value_t>(matrix(i, j)) << " ";
      }
      std::cout << std::endl;
    }
  }
}


template<typename MatrixT>
static
void verify_matrix(MatrixT &result, MatrixT &expected)
{
  using LocalBlockCache = MatrixBlock<MatrixT>;
  const size_t num_blocks = result.pattern().blockspec().extent(0);
  for (size_t i = 0; i < num_blocks; ++i) {
    for (size_t j = 0; j < num_blocks; ++j) {
      LocalBlockCache block_r(result, i, j);
      if (block_r.is_local()) {
        LocalBlockCache block_e(expected, i, j);
        auto e = block_e.lbegin();
        auto r = block_r.lbegin();
        for (size_t k = 0; k < block_r.size(); ++k) {
          if (e[k] != r[k]) {
            std::cout << "Result diverges in block " << i << "x" << j << " by "
                      << e[k] - r[k] << "(" << e[k] << " vs " << r[k] << ")" << std::endl;
            //return;
          }
        }
      }
    }
  }
  std::cout << "Verification successful!" << std::endl;
}



#endif // CHOLESKY_COMMON_H
