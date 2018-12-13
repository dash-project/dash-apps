#ifndef CHOLESKY_TASKS_PREFETCH__H
#define CHOLESKY_TASKS_PREFETCH__H

#include <libdash.h>
#include <memory>
#include <thread>
#include "MatrixBlock.h"
#include "common.h"

constexpr const char *QRFACT_IMPL = "QRFactTasksCopyin";

/**
 * TODO: this implementation only supports quadratic tiles because MatrixBlock
 *       is a very simple wrapper. Fix this as soon as DASH matrix blocks
 *       are actually usable.
 */

thread_local double *scratch_tau  = nullptr;
thread_local double *scratch_work = nullptr;

template<typename TiledMatrix>
void
compute(TiledMatrix& A, TiledMatrix& T, size_t block_size){

  using value_t       = typename TiledMatrix::value_type;
  using Block         = MatrixBlock<TiledMatrix>;
  using BlockCache    = typename std::vector<value_t>;
  using BlockCachePtr = typename std::shared_ptr<BlockCache>;
  const size_t num_blocks = A.pattern().extent(0) / block_size;

  // allocate a vector to pre-fetch the result of trsm() into
  value_t *a_block_k_pre   = new value_t[block_size*block_size];
  value_t *t_block_k_pre   = new value_t[block_size*block_size];
  value_t *a_blocks_mk_pre = new value_t[block_size*block_size*num_blocks];
  value_t *t_blocks_mk_pre = new value_t[block_size*block_size*num_blocks];
  value_t *a_blocks_kn_pre = new value_t[block_size*block_size*num_blocks];

  auto next_buf_pos = [=](){
    static size_t buffer_pos = 0;
    size_t res = buffer_pos;
    buffer_pos = (buffer_pos + 1) % num_blocks;
    return res*block_size*block_size;
  };

  size_t num_tasks = 0;
  if (dash::myid() == 0)
    std::cout << "Starting to create tasks" << std::endl;

  Timer t_c;

  // iterate over column of blocks
  for (int k = 0; k < num_blocks; ++k) {

    Block a_block_k(A, k, k);
    Block t_block_k(T, k, k);
    value_t *a_k_pre;
    value_t *t_k_pre;

    /**
     * Compute cholesky factorization of block on diagonal
     */
    if (a_block_k.is_local()) {
      dash::tasks::async("GEQRT",
        [=]() mutable {
          if (scratch_work == nullptr) scratch_work = new double[block_size*block_size];
          if (scratch_tau  == nullptr) scratch_tau  = new double[block_size];
          geqrt(a_block_k.lbegin(), t_block_k.lbegin(), scratch_tau, scratch_work, block_size);
#ifdef DEBUG_OUTPUT
          printf("\n\nPOST geqrt(k=%d)\n", k);
          print_matrix(a_block_k.lbegin(), block_size);
          print_matrix(t_block_k.lbegin(), block_size);
#endif
        },
        //(dart_task_prio_t)((num_blocks-k)*(num_blocks-k)*(num_blocks-k))/*priority*/,
        dash::tasks::out(a_block_k),
        dash::tasks::out(t_block_k)
      );
      ++num_tasks;

      a_k_pre = a_block_k.lbegin();
      t_k_pre = t_block_k.lbegin();
    } else {
      a_k_pre = a_block_k_pre;
      t_k_pre = t_block_k_pre;
    }
    dash::tasks::async_barrier();

    for (int n = k+1; n < num_blocks; ++n) {
      Block a_block_kn(A, k, n);
      if (a_block_kn.is_local()) {
        dash::tasks::async("ORMQR",
          [=]() mutable {
            if (scratch_work == nullptr) scratch_work = new double[block_size*block_size];
#ifdef DEBUG_OUTPUT
            printf("\n\nPRE ormqr(k=%d, n=%d)\n", k, n);
            print_matrix(a_k_pre, block_size);
            print_matrix(t_k_pre, block_size);
            print_matrix(a_block_kn.lbegin(), block_size);
#endif
            ormqr(a_k_pre, t_k_pre, a_block_kn.lbegin(), scratch_work, block_size);

#ifdef DEBUG_OUTPUT
            printf("\n\nPOST ormqr(k=%d, n=%d)\n", k, n);
            print_matrix(a_block_kn.lbegin(), block_size);
#endif
          },
          DART_PRIO_LOW,
          dash::tasks::copyin_r(a_block_k, block_size*block_size, a_k_pre),
          dash::tasks::copyin_r(t_block_k, block_size*block_size, t_k_pre),
          dash::tasks::out(a_block_kn)
        );
        ++num_tasks;
      }
    }

    //dash::tasks::async_barrier();

    for (int m = k+1; m < num_blocks; ++m) {

      Block a_block_mk(A, m, k);
      Block t_block_mk(T, m, k);
      if (a_block_mk.is_local()) {
        dash::tasks::async("TSQRT",
          [=]() mutable {
            if (scratch_work == nullptr) scratch_work = new double[block_size*block_size];
            if (scratch_tau  == nullptr) scratch_tau  = new double[block_size];
            tsqrt(a_k_pre,
                  a_block_mk.lbegin(),
                  t_block_mk.lbegin(),
                  scratch_tau, scratch_work, block_size);
#ifdef DEBUG_OUTPUT
            printf("\n\ntsqrt(k=%d, m=%d)\n", k, m);
            print_matrix(a_k_pre, block_size);
            print_matrix(a_block_mk.lbegin(), block_size);
            print_matrix(t_block_mk.lbegin(), block_size);
#endif

            if (!a_block_k.is_local()) {
              a_block_k.store_async(a_k_pre);
            }
            while (!a_block_k.test()) dash::tasks::yield(-1);
          },
          //(dart_task_prio_t)((num_blocks-k)*(num_blocks-k)*(num_blocks-k)) /*priority*/,
          dash::tasks::copyin_r(a_block_k, block_size*block_size, a_k_pre),
          dash::tasks::out(a_block_mk),
          dash::tasks::out(t_block_mk),
          dash::tasks::out(a_block_k)
        );
        ++num_tasks;
      }

      dash::tasks::async_barrier();

      value_t *a_block_mk_pre;
      value_t *t_block_mk_pre;
      if (a_block_mk.is_local()) {
        a_block_mk_pre = a_block_mk.lbegin();
        t_block_mk_pre = t_block_mk.lbegin();
      } else {
        a_block_mk_pre = &a_blocks_mk_pre[block_size*block_size*m];
        t_block_mk_pre = &t_blocks_mk_pre[block_size*block_size*m];
      }

      for (int n = k+1; n < num_blocks; ++n) {

        Block a_block_mn(A, m, n);
        Block a_block_kn(A, k, n);
        value_t *a_block_kn_pre;
        if (a_block_kn.is_local()) {
          a_block_kn_pre = a_block_kn.lbegin();
        } else {
          a_block_kn_pre = &a_blocks_kn_pre[block_size*block_size*m];
        }

        if (a_block_mn.is_local()) {
          dash::tasks::async("TSMQR",
            [=]() mutable {
              if (scratch_work == nullptr) scratch_work = new double[block_size*block_size];

              tsmqr(a_block_kn_pre,
                    a_block_mn.lbegin(),
                    a_block_mk_pre,
                    t_block_mk_pre,
                    scratch_work, block_size);

#ifdef DEBUG_OUTPUT
              printf("\n\ntsmqr(k=%d, m=%d)\n", k, m);
              print_matrix(a_block_kn_pre, block_size);
              print_matrix(a_block_mn.lbegin(), block_size);
              print_matrix(a_block_mk_pre, block_size);
              print_matrix(t_block_mk_pre, block_size);
#endif

              // store the block that might not be local
              if (!a_block_kn.is_local()) {
                a_block_kn.store_async(a_block_kn_pre);
              }
              while (!a_block_kn.test()) dash::tasks::yield(-1);
            },
            //(dart_task_prio_t)((num_blocks-k)*(num_blocks-n)*(num_blocks-n)),
            dash::tasks::copyin_r(a_block_kn, block_size*block_size, a_block_kn_pre),
            dash::tasks::copyin_r(a_block_mk, block_size*block_size, a_block_mk_pre),
            dash::tasks::copyin_r(t_block_mk, block_size*block_size, t_block_mk_pre),
            dash::tasks::out(a_block_mn),
            dash::tasks::out(a_block_kn)
          );
          ++num_tasks;
        }
      }
    }
    dash::tasks::async_barrier();
  }

  if (dash::myid() == 0)
    std::cout << "Done creating " << num_tasks << " tasks in "
              << t_c.Elapsed() / 1E3 << "ms"
              << ", starting execution" << std::endl;
  dash::tasks::complete();
  if (dash::myid() == 0)
    std::cout << "Done executing " << num_tasks << " tasks after "
              << t_c.Elapsed() / 1E3 << "ms"
              << std::endl;

  delete[] a_block_k_pre;
  delete[] t_block_k_pre;
  delete[] a_blocks_mk_pre;
  delete[] t_blocks_mk_pre;
  delete[] a_blocks_kn_pre;

}

#endif //CHOLESKY_TASKS_PREFETCH__H
