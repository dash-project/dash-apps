#ifndef CHOLESKY_COMMON_H
#define CHOLESKY_COMMON_H

// TODO: core_cblas.h is not self-contained?
#include <plasma.h>
#include <core_dblas.h>
#include "MatrixBlock.h"
#include "ExtraeInstrumentation.h"

typedef dash::util::Timer<dash::util::TimeMeasure::Clock> Timer;

template<typename ValueT>
static
void print_matrix(ValueT *matrix, size_t num_elem)
{
  if (num_elem > 100) return;

  using value_t = ValueT;

  std::cout << std::setw(8);

  std::cout << "###### Matrix(" << num_elem << "," << num_elem << ") ######" << std::endl;

  if (dash::myid() == 0) {
    for (size_t i = 0; i < num_elem; ++i) {
      for (size_t j = 0; j < num_elem; ++j) {
        std::cout << std::setw(8) << std::setprecision(4) << std::fixed
                  << static_cast<value_t>(matrix[i*num_elem + j]) << " ";
      }
      std::cout << std::endl;
    }
  }
  std::cout << "###############" << std::endl;
}

template<typename PatternT, typename ValueT>
static
void print_matrix(dash::Matrix<ValueT, 4, dash::default_index_t, PatternT>& matrix)
{
  using value_t = ValueT;

  std::cout << std::setw(8);

  std::cout << "!###### Matrix(" << matrix.pattern().extent(0) << ") ######" << std::endl;

  if (dash::myid() == 0) {
    for (size_t i = 0; i < matrix.pattern().extent(0); ++i) {
      for (size_t k = 0; k < matrix.pattern().extent(2); ++k) {
        for (size_t j = 0; j < matrix.pattern().extent(1); ++j) {
          for (size_t l = 0; l < matrix.pattern().extent(3); ++l) {
            std::cout << std::setw(8) << std::setprecision(4) << std::fixed
                << static_cast<ValueT>(matrix(i, j, k, l)) << " ";
          }
        }
        std::cout << std::endl;
      }
    }
  }

  std::cout << "###############" << std::endl;
}



template<typename ValueT>
static
void print_block(ValueT *block, size_t num_elem)
{
  if (num_elem > 100) return;

  using value_t = ValueT;

  std::cout << std::setw(8);

  std::cout << "###### Matrix(" << num_elem << "," << num_elem << ") " << block << " ######" << std::endl;

  for (size_t i = 0; i < num_elem; ++i) {
    for (size_t j = 0; j < num_elem; ++j) {
      std::cout << std::setw(8) << std::setprecision(4) << std::fixed
                << static_cast<value_t>(block[i*num_elem + j]) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "###############" << std::endl;
}

#if 0

template<typename MatrixT>
static
void print_matrix(MatrixT &matrix)
{
  if (matrix.extent(0) > 100) return;

  using value_t = typename MatrixT::value_type;

  std::cout << std::setw(8);

  if (dash::myid() == 0) {
    // 4-dimensional matrix modeling a two-dimensional matrix with super-blocks
    size_t n_row = matrix.extent(0) * matrix.extent(2);
    size_t c = 0;
    for (auto it = matrix.begin(); it != matrix.end(); ++it) {
      std::cout << std::setw(8) << std::setprecision(4) << std::fixed
                << static_cast<value_t>(*it) << " ";
      if (++c % n_row == 0) std::cout << std::endl;
    }
  }
}

#endif

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


static
void geqrt(double *A, double *T, double *TAU, double *WORK, int ts)
{
  //std::cout << "BEFORE dgeqrt: A:\n";
  //print_block(A, ts);
  //std::cout << "BEFORE dgeqrt: T:\n";
  //print_block(T, ts);
  int ret = CORE_dgeqrt(ts, ts, ts, A, ts, T, ts, TAU, WORK);
  assert(ret == 0);
  //std::cout << "dgeqrt: A:\n";
  //print_block(A, ts);
  //std::cout << "dgeqrt: T:\n";
  //print_block(T, ts);
}

static
void ormqr(double *A, double *T, double *C, double *WORK, int ts)
{
  int side = PlasmaLeft;
  int trans = PlasmaTrans;

  //std::cout << "BEFORE dormqr: A:\n";
  //print_block(A, ts);
  //std::cout << "BEFORE dormqr: T:\n";
  //print_block(T, ts);
  //std::cout << "BEFORE dormqr: C:\n";
  //print_block(C, ts);
  int ret = CORE_dormqr(side, trans, ts, ts, ts, ts,
                  A, ts, T, ts, C, ts, WORK, ts);
  assert(ret == 0);
  //std::cout << "dormqr: A:\n";
  //print_block(A, ts);
  //std::cout << "dormqr: T:\n";
  //print_block(T, ts);
  //std::cout << "dormqr: C:\n";
  //print_block(C, ts);
}

static
void tsqrt(double *A, double *B, double *T, double *TAU, double *WORK, int ts)
{
  //std::cout << "BEFORE dtsqrt: A:\n";
  //print_block(A, ts);

  //std::cout << "BEFORE dtsqrt: B:\n";
  //print_block(B, ts);

  //std::cout << "BEFORE dtsqrt: T:\n";
  //print_block(T, ts);

  int ret = CORE_dtsqrt(ts, ts, ts, A, ts, B, ts, T, ts, TAU, WORK);
  assert(ret == 0);
  //std::cout << "dtsqrt: A:\n";
  //print_block(A, ts);

}

static
void tsmqr(double *A, double *B, double *V, double *T, double* WORK, int ts)
{
  int side = PlasmaLeft;
  int trans = PlasmaTrans;

  //std::cout << "BEFORE dtsmqr: A:\n";
  //print_block(A, ts);
  int ret = CORE_dtsmqr(side, trans, ts, ts, ts, ts, ts, ts,
            A, ts, B, ts, V, ts, T, ts, WORK, ts);
  assert(ret == 0);
  //std::cout << "dtsmqr: A:\n";
  //print_block(A, ts);

}

template<typename MatrixT>
static
void
compute_single(MatrixT& A, MatrixT& T, size_t block_size){

  using LocalBlockCache = MatrixBlock<MatrixT>;


  using value_t       = typename MatrixT::value_type;
  const size_t num_blocks = A.pattern().extent(0) / block_size;

  double *work = new double[block_size*block_size];
  double *tau  = new double[block_size];

  // iterate over column of blocks
  for (int k = 0; k < num_blocks; ++k) {

    LocalBlockCache a_block_k(A, k, k);
    LocalBlockCache t_block_k(T, k, k);

    if (a_block_k.is_local()) {
      printf("\n\ngeqrt(k=%d)\n", k);
      print_matrix(a_block_k.lbegin(), block_size);
      print_matrix(t_block_k.lbegin(), block_size);
      geqrt(a_block_k.lbegin(), t_block_k.lbegin(), tau, work, block_size);
      printf("\n\ngeqrt(k=%d)\n", k);
      print_matrix(a_block_k.lbegin(), block_size);
      print_matrix(t_block_k.lbegin(), block_size);
    }
    dash::barrier();

    if (!a_block_k.is_local()) {
      a_block_k.fetch();
      t_block_k.fetch();
    }

    for (int n = k+1; n < num_blocks; ++n) {
      LocalBlockCache a_block_kn(A, k, n);
      if (a_block_kn.is_local()) {
        ormqr(a_block_k.lbegin(), t_block_k.lbegin(),
              a_block_kn.lbegin(), work, block_size);
        printf("\n\normqr(k=%d, n=%d)\n", k, n);
        print_matrix(a_block_kn.lbegin(), block_size);
      }
    }
    dash::barrier();

    for (int m = k+1; m < num_blocks; ++m) {

      LocalBlockCache a_block_mk(A, m, k);
      LocalBlockCache t_block_mk(T, m, k);
      if (a_block_mk.is_local()) {
        if (!a_block_k.is_local()) {
          a_block_k.fetch();
        }

        tsqrt(a_block_k.lbegin(),
              a_block_mk.lbegin(),
              t_block_mk.lbegin(),
              tau, work, block_size);

        printf("\n\ntsqrt(k=%d, m=%d)\n", k, m);
        print_matrix(a_block_k.lbegin(), block_size);
        print_matrix(a_block_mk.lbegin(), block_size);
        print_matrix(t_block_mk.lbegin(), block_size);

        if (!a_block_k.is_local()) {
          a_block_k.store();
        }
      }

      dash::barrier();

      if (!a_block_mk.is_local()) {
        a_block_mk.fetch();
        t_block_mk.fetch();
      }

      for (int n = k+1; n < num_blocks; ++n) {

        LocalBlockCache a_block_mn(A, m, n);
        LocalBlockCache a_block_kn(A, k, n);

        if (a_block_mn.is_local()) {
          if (!a_block_kn.is_local()) {
            a_block_kn.fetch();
          }

          tsmqr(a_block_kn.lbegin(),
                a_block_mn.lbegin(),
                a_block_mk.lbegin(),
                t_block_mk.lbegin(),
                work, block_size);

          printf("\n\ntsmqr(k=%d, m=%d)\n", k, m);
          print_matrix(a_block_kn.lbegin(), block_size);
          print_matrix(a_block_mn.lbegin(), block_size);
          print_matrix(a_block_mk.lbegin(), block_size);
          print_matrix(t_block_mk.lbegin(), block_size);

          // store the block that might not be local
          if (!a_block_kn.is_local()) {
            a_block_kn.store();
          }
        }
      }
      dash::barrier();
    }
  }

  delete[] tau;
  delete[] work;
}

#endif // CHOLESKY_COMMON_H
