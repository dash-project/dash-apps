#ifndef CHOLESKY_TASKS_PREFETCH__H
#define CHOLESKY_TASKS_PREFETCH__H

#include <libdash.h>
#include <memory>
#include "MatrixBlock.h"
#include "common.h"
#include "BlockPrefetcher.h"

constexpr const char *CHOLESKY_IMPL = "CholeskyTasksCopyin";

//#define DEBUG_OUT

template<typename TiledMatrix>
void
compute(TiledMatrix& matrix, int block_size){

  using value_t       = typename TiledMatrix::value_type;
  using Block         = MatrixBlock<TiledMatrix>;
  using BlockCache    = typename std::vector<value_t>;
  using BlockCachePtr = typename std::shared_ptr<BlockCache>;
  const size_t num_blocks = matrix.pattern().extent(0) / block_size;

  // allocate a vector to pre-fetch the result of trsm() into
  value_t *blocks_ki_pre = new value_t[block_size*block_size*num_blocks];

  /**
   * Algorithm taken from
   * https://pm.bsc.es/ompss-docs/examples/01-examples/cholesky
   */

  size_t num_tasks = 0;
  if (dash::myid() == 0)
    std::cout << "Starting to create tasks" << std::endl;

  Timer t_c;

  int buffer_pos = 0;
  // iterate over column of blocks
  for (int k = 0; k < num_blocks; ++k) {

    Block block_k(matrix, k, k);

#ifdef DEBUG_OUT
    if (dash::myid() == 0)
      std::cout << "k = " << k << std::endl;
#endif

    /**
     * Compute cholesky factorization of block on diagonal
     */
    if (block_k.is_local()) {
      dash::tasks::async("POTRF",
        [&matrix, k, block_size, num_blocks/*, block_k*/]() {
#ifdef DEBUG_OUT
          std::cout << "[" << dash::myid() << ", " << dart_task_thread_num()
                    << "] potrf() on row " << k << "/" << num_blocks << ": ";
#endif
          Block block_k(matrix, k, k);
          potrf(block_k.lbegin(), block_size, block_size);
#ifdef DEBUG_OUT
          std::cout << "Done." << std::endl;
#endif
        },
        //DART_PRIO_HIGH,
        (dart_task_prio_t)((num_blocks - k) * (num_blocks-k) * (num_blocks - k))/*priority*/,
        dash::tasks::out(block_k)
      );
      ++num_tasks;
    }
    dash::tasks::async_fence();

    /**
     * Solve the triangular equation system in the block
     * Only block_k needs prefetching here.
     */
    for (int i = k+1; i < num_blocks; ++i) {
      Block block_b(matrix, k, i);
      if (block_b.is_local()) {
        //auto block_k_pre = &blocks_ki_pre[k*block_size*block_size];
        auto block_k_pre = block_k.is_local()
                              ? block_k.lbegin() : &blocks_ki_pre[buffer_pos*block_size*block_size];
        dash::tasks::async("TRSM",
          [&matrix, block_k_pre, k, i, block_size/*, block_b*/]() {
            Block block_b(matrix, k, i);
            trsm(block_k_pre,
                block_b.lbegin(), block_size, block_size);
          },
          //DART_PRIO_HIGH,
          (dart_task_prio_t)((num_blocks - i) * (num_blocks-i) * (num_blocks - i) + 3 * ((2 * num_blocks) - k - i - 1) * (i - k))/*priority*/,
          dash::tasks::copyin_r(block_k.begin(), block_k.end(), block_k_pre),
          dash::tasks::out(block_b)
        );
        ++num_tasks;
      }
    }
    dash::tasks::async_fence();

    buffer_pos = (buffer_pos+1)%num_blocks;

    // walk to the right
    int ib = buffer_pos;
    for (int i = k+1; i < num_blocks; ++i, ib = (ib+1)%num_blocks) {
      Block block_a(matrix, k, i);
      int i_buffer_pos = (buffer_pos+1) % num_blocks;
      auto block_a_pre = block_a.is_local()
                            ? block_a.lbegin()
                            : &blocks_ki_pre[ib*block_size*block_size];

      // update diagonal blocks
      Block block_i(matrix, i, i);
      if (block_i.is_local()) {
        // A[j,j] = A[j,j] - A[j,i] * (A[j,i])^t
        dash::tasks::async("SYRK",
          [&matrix, i, block_size, block_a_pre/*, block_i*/]() mutable {
            Block block_i(matrix, i, i);
            syrk(block_a_pre,
                 block_i.lbegin(), block_size, block_size);
          },
          //DART_PRIO_HIGH,
          (dart_task_prio_t)((num_blocks - i) * (num_blocks - i) * (num_blocks - i) + 3 * (i - k))/*priority*/,
          dash::tasks::copyin_r(block_a.begin(), block_a.end(), block_a_pre),
          dash::tasks::out(block_i)
        );
        ++num_tasks;
      }

      // run down to the diagonal
      int jb = buffer_pos;
      for (int j = i+1; j < num_blocks; ++j, jb = (jb+1)%num_blocks) {
        Block block_c(matrix, i, j);
        if (block_c.is_local()) {
          Block block_b(matrix, k, j);
          auto block_b_pre = block_b.is_local()
                              ? block_b.lbegin()
                              : &blocks_ki_pre[jb*block_size*block_size];
          // A[k,i] = A[k,i] - A[k,j] * (A[j,i])^t
          dash::tasks::async("GEMM",
            [&matrix, block_a_pre, block_b_pre, block_size, i, j/*, block_c*/]() mutable {
              Block block_c(matrix, i, j);
              gemm(block_a_pre, block_b_pre,
                   block_c.lbegin(), block_size, block_size);
            },
            (dart_task_prio_t)((num_blocks - i) * (num_blocks - i) * (num_blocks - i) + 3 * ((2 * num_blocks) - i - j - 3) * (i - j) + 6 * (i - k)) /*priority*/,
            dash::tasks::copyin_r(block_a.begin(), block_a.end(), block_a_pre),
            dash::tasks::copyin_r(block_b.begin(), block_b.end(), block_b_pre),
            dash::tasks::out(block_c)
          );
          ++num_tasks;
        }
      }
    }
    dash::tasks::async_fence();
    // ib points to next block not used by loop above
    buffer_pos = ib;
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

  delete[] blocks_ki_pre;

}

#endif //CHOLESKY_TASKS_PREFETCH__H
