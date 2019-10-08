#ifndef CHOLESKY_TASKS_PREFETCH__H
#define CHOLESKY_TASKS_PREFETCH__H

#include <libdash.h>
#include <memory>
#include <thread>
#include <vector>
#include "MatrixBlock.h"
#include "common.h"

constexpr const char *QRFACT_IMPL = "QRFactTasksCopyin";

/**
 * TODO: this implementation only supports quadratic tiles because MatrixBlock
 *       is a very simple wrapper. Fix this as soon as DASH matrix blocks
 *       are actually usable.
 */

template<typename TiledMatrix>
void
compute(TiledMatrix& A, TiledMatrix& T, size_t block_size){

  using value_t       = typename TiledMatrix::value_type;
  using Block         = MatrixBlock<TiledMatrix>;
  using BlockCache    = typename std::vector<value_t>;
  using BlockCachePtr = typename std::shared_ptr<BlockCache>;
  const size_t num_blocks = A.pattern().extent(0);

  // allocate a vector to pre-fetch the result of trsm() into
  value_t *a_block_k_pre   = new value_t[block_size*block_size];
  value_t *t_block_k_pre   = new value_t[block_size*block_size];
  value_t *a_blocks_mk_pre = new value_t[block_size*block_size*num_blocks];
  value_t *t_blocks_mk_pre = new value_t[block_size*block_size*num_blocks];

  dash::Matrix<value_t, 3> scratch_blocks_kn(dash::size(), dash::BLOCKED, num_blocks, dash::NONE, block_size*block_size, dash::NONE);
  constexpr const int num_kk_scratch = 2;
  dash::Matrix<value_t, 3> scratch_blocks_kk(dash::size(), dash::BLOCKED, num_kk_scratch, dash::NONE, block_size*block_size, dash::NONE);

  double **scratch_tau  = new double*[dash::tasks::numthreads()];
  std::fill(scratch_tau, scratch_tau+dash::tasks::numthreads(), nullptr);
  double **scratch_work = new double*[dash::tasks::numthreads()];
  std::fill(scratch_work, scratch_work+dash::tasks::numthreads(), nullptr);

  size_t num_tasks = 0;
  if (dash::myid() == 0)
    std::cout << "Starting to create tasks" << std::endl;

  Timer t_c;

  std::cout << "A: " << A.lbegin() << " - " << A.lend() << std::endl;

#if 0
  for (int i = 0; i < num_blocks; i++) {
    for (int j = 0; j < num_blocks; j++) {
      std::cout << "Block {" << i << ", " << j << "}" << std::endl;
      print_block(Block{A, i, j}.lbegin(), block_size);
    }
  }
#endif

  // iterate over column of blocks
  for (int k = 0; k < num_blocks; ++k) {

    Block a_block_k(A, k, k);
    Block t_block_k(T, k, k);
    value_t *a_k_pre;
    value_t *t_k_pre;
    int k_scratch_entry = k%num_kk_scratch;

    /**
     * Compute cholesky factorization of block on diagonal
     */
    if (a_block_k.is_local()) {
      dash::tasks::async("GEQRT",
        [=, &scratch_blocks_kk, &A, &T, &scratch_work, &scratch_tau]() mutable {
          int threadnum = dash::tasks::threadnum();
          if (scratch_work[threadnum] == nullptr) scratch_work[threadnum] = new double[block_size*block_size];
          if (scratch_tau[threadnum]  == nullptr) scratch_tau[threadnum]  = new double[block_size];
          Block a_block_k(A, k, k);
          Block t_block_k(T, k, k);
          std::cout << "geqrt " << k << ", " << k << std::endl;
          geqrt(a_block_k.lbegin(), t_block_k.lbegin(), scratch_tau[threadnum], scratch_work[threadnum], block_size);

          std::copy(a_block_k.lbegin(), a_block_k.lend(),
                    &scratch_blocks_kk.local(0, k_scratch_entry, 0));
#ifdef DEBUG_OUTPUT
          printf("\n\nPOST geqrt(k=%d)\n", k);
          print_matrix(a_block_k.lbegin(), block_size);
          print_matrix(t_block_k.lbegin(), block_size);
#endif
        },
        (dart_task_prio_t)((num_blocks-k)*(num_blocks-k)*(num_blocks-k))/*priority*/,
        dash::tasks::out(a_block_k),
        dash::tasks::out(t_block_k),
        dash::tasks::out(scratch_blocks_kk.local(0, k_scratch_entry, 0))
      );
      ++num_tasks;

      a_k_pre = a_block_k.lbegin();
      t_k_pre = t_block_k.lbegin();
    } else {
      a_k_pre = a_block_k_pre;
      t_k_pre = t_block_k_pre;
    }
    dash::tasks::async_fence();

    for (int n = k+1; n < num_blocks; ++n) {
      Block a_block_kn(A, k, n);
      if (a_block_kn.is_local()) {
        dash::tasks::async("ORMQR",
          [=, &scratch_blocks_kn, &A, &scratch_work]() mutable {
            int threadnum = dash::tasks::threadnum();
            if (scratch_work[threadnum] == nullptr) scratch_work[threadnum] = new double[block_size*block_size];
#ifdef DEBUG_OUTPUT
            printf("\n\nPRE ormqr(k=%d, n=%d)\n", k, n);
            print_matrix(a_k_pre, block_size);
            print_matrix(t_k_pre, block_size);
            print_matrix(a_block_kn.lbegin(), block_size);
            std::cout << "ormqr " << k << ", " << n << std::endl;
            std::cout << a_block_kn.lbegin() << std::endl;
#endif
            Block a_block_kn(A, k, n);
            // Block(k, n) is both input and output, use the original block as
            // input and copy it to the scratch space afterwards
            ormqr(a_k_pre, t_k_pre, a_block_kn.lbegin(), scratch_work[threadnum], block_size);
            // copy the block to the scratch space for subsequent accesses
            std::copy(a_block_kn.lbegin(), a_block_kn.lend(),
                      &scratch_blocks_kn.local(0, n, 0));
#ifdef DEBUG_OUTPUT
            printf("\n\nPOST ormqr(k=%d, n=%d)\n", k, n);
            print_matrix(a_block_kn.lbegin(), block_size);
#endif
          },
          DART_PRIO_LOW,
          dash::tasks::copyin_r(a_block_k, block_size*block_size, a_k_pre),
          dash::tasks::copyin_r(t_block_k, block_size*block_size, t_k_pre),
          dash::tasks::out(a_block_kn),
          dash::tasks::out(scratch_blocks_kn.local(0, n, 0))
        );
        ++num_tasks;
      }
    }

    //dash::tasks::async_fence();

    for (int m = k+1; m < num_blocks; ++m) {

      Block a_block_mk(A, m, k);
      Block t_block_mk(T, m, k);
      if (a_block_mk.is_local()) {
        value_t *a_block_kk_pre;
        // there was a previous owner, get it directly from his scratch space
        a_block_kk_pre = &scratch_blocks_kk.local(0, k_scratch_entry, 0);
        int prev_k_owner = Block{A, m-1, k}.unit();
        dash::tasks::async("TSQRT",
          [=, &A, &T, &scratch_work, &scratch_tau]() mutable {
            int threadnum = dash::tasks::threadnum();
            if (scratch_work[threadnum] == nullptr) scratch_work[threadnum] = new double[block_size*block_size];
            if (scratch_tau[threadnum]  == nullptr) scratch_tau[threadnum]  = new double[block_size];
#ifdef DEBUG_OUTPUT
            printf("\n\nPRE tsqrt(k=%d, n=%d)\n", k, n);
            print_matrix(a_k_pre, block_size);
            print_matrix(t_k_pre, block_size);
            print_matrix(a_block_kn.lbegin(), block_size);
            std::cout << "tsqrt " << k << ", " << m << std::endl;
#endif
            Block a_block_mk(A, m, k);
            Block t_block_mk(T, m, k);
            tsqrt(a_block_kk_pre,
                  a_block_mk.lbegin(),
                  t_block_mk.lbegin(),
                  scratch_tau[threadnum], scratch_work[threadnum], block_size);
#ifdef DEBUG_OUTPUT
            printf("\n\nPOST tsqrt(k=%d, m=%d)\n", k, m);
            print_matrix(a_block_kk_pre, block_size);
            print_matrix(a_block_mk.lbegin(), block_size);
            print_matrix(t_block_mk.lbegin(), block_size);
#endif

          },
          (dart_task_prio_t)((num_blocks-k)*(num_blocks-k)*(num_blocks-k)) /*priority*/,
          dash::tasks::copyin_r(scratch_blocks_kk(prev_k_owner, k_scratch_entry, 0), block_size*block_size, a_block_kk_pre),
          dash::tasks::out(a_block_mk),
          dash::tasks::out(t_block_mk),
          dash::tasks::out(scratch_blocks_kk.local(0, k_scratch_entry, 0))
        );
        ++num_tasks;
      }

      dash::tasks::async_fence();

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
        if (a_block_mn.is_local()) {
          value_t *a_block_kn_pre;
          // there was a previous owner, get it directly from his scratch space
          a_block_kn_pre = &scratch_blocks_kn.local(0, n, 0);
          int prev_owner = Block{A, m-1, n}.unit();
          dash::tasks::async("TSMQR",
            [=, &A, &scratch_work]() mutable {
              int threadnum = dash::tasks::threadnum();
              if (scratch_work[threadnum] == nullptr) scratch_work[threadnum] = new double[block_size*block_size];
#ifdef DEBUG_OUTPUT
              std::cout << "tsmqr " << k << ", " << m << ", " << n << std::endl;
#endif
              Block a_block_mn(A, m, n);
              tsmqr(a_block_kn_pre,
                    a_block_mn.lbegin(),
                    a_block_mk_pre,
                    t_block_mk_pre,
                    scratch_work[threadnum], block_size);

#ifdef DEBUG_OUTPUT
              printf("\n\ntsmqr(k=%d, m=%d)\n", k, m);
              print_matrix(a_block_kn_pre, block_size);
              print_matrix(a_block_mn.lbegin(), block_size);
              print_matrix(a_block_mk_pre, block_size);
              print_matrix(t_block_mk_pre, block_size);
#endif
            },
            (dart_task_prio_t)((num_blocks-k)*(num_blocks-n)*(num_blocks-n)),
            dash::tasks::copyin_r(scratch_blocks_kn(prev_owner, n, 0), block_size*block_size, a_block_kn_pre),
            dash::tasks::copyin_r(a_block_mk, block_size*block_size, a_block_mk_pre),
            dash::tasks::copyin_r(t_block_mk, block_size*block_size, t_block_mk_pre),
            dash::tasks::out(a_block_mn)
          );
          ++num_tasks;
        }
      }
    }
    dash::tasks::async_fence();

    // write back block (k,k)
    int prev_k_owner = Block{A, num_blocks-1, k}.unit();
    if (scratch_blocks_kk(prev_k_owner, k_scratch_entry, 0).is_local()) {
      // write back Block (k, k)
      dash::tasks::async("WRITEBACK_KK",
        [=, &scratch_blocks_kk, &A]() mutable {
            // write the block back
            Block a_block_k(A, k, k);
            a_block_k.store_async(&scratch_blocks_kk.local(0, k_scratch_entry, 0));
            // detach the task but release the dependencies only when the transfer completed
            dart_handle_t handle = a_block_k.dart_handle();
            dart_task_detach_handle(&handle, 1);
        },
        (dart_task_prio_t)((num_blocks-k)*(num_blocks-k)*(num_blocks-k)) /*priority*/,
        dash::tasks::in(scratch_blocks_kk.local(0, k_scratch_entry, 0))
        // NOTE: we don't need an output dependency here as block (k, k) will never be touched again'
        //dash::tasks::out(a_block_k)
      );
    }

    // write back blocks (k, n)
    for (int n = k+1; n < num_blocks; ++n) {
      int prev_owner = Block{A, num_blocks-1, n}.unit();
      auto scratch_block_kn = scratch_blocks_kn(prev_owner, n, 0);
      if (scratch_block_kn.is_local()) {
        dash::tasks::async("WRITEBACK_KN",
          [=, &scratch_blocks_kn, &A]() mutable {
            Block a_block_kn(A, k, n);
            // write the block back
            a_block_kn.store_async(&scratch_blocks_kn.local(0, n, 0));
            // detach the task but release the dependencies only when the transfer completed
            dart_handle_t handle = a_block_kn.dart_handle();
            dart_task_detach_handle(&handle, 1);
          },
          (dart_task_prio_t)((num_blocks-k)*(num_blocks-n)*(num_blocks-n)),
          dash::tasks::in(scratch_block_kn)
          // NOTE: we don't need an output dependency because block(k, n)'is never again read
        );
        ++num_tasks;
      }
    }
  }

  if (dash::myid() == 0)
    std::cout << "Done creating " << num_tasks << " tasks in "
              << t_c.Elapsed() / 1E3 << "ms"
              << ", starting execution" << std::endl;
  dash::tasks::complete();

  auto elapsed = t_c.Elapsed();
  size_t num_elem = A.extent(0);
  double flops = (4.0 / 3.0) * num_elem * num_elem * num_elem;
  if (dash::myid() == 0)
    std::cout << "Done executing " << num_tasks << " tasks after "
              << t_c.Elapsed() / 1E3 << "ms"
              << std::endl;

/*
  if (dash::myid() == 0) {
    std::cout << "QR factorization of " << num_elem << "x" << num_elem
              << " done after " << elapsed/1E3 << "ms ("
              << flops/elapsed/1E3 << " GF/s)" << std::endl;
  }
*/

  delete[] a_block_k_pre;
  delete[] t_block_k_pre;
  delete[] a_blocks_mk_pre;
  delete[] t_blocks_mk_pre;

  for (int i = 0; i < dash::tasks::numthreads(); ++i) {
    if (scratch_tau[i]  != nullptr) delete[] scratch_tau[i];
    if (scratch_work[i] != nullptr) delete[] scratch_work[i];
  }
  delete[] scratch_tau;
  delete[] scratch_work;
}

#endif //CHOLESKY_TASKS_PREFETCH__H
