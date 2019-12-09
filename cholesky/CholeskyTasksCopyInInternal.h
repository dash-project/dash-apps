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

  /**
   * Algorithm taken from
   * https://pm.bsc.es/ompss-docs/examples/01-examples/cholesky
   */

  size_t num_tasks = 0;
  if (dash::myid() == 0)
    std::cout << "Starting to create tasks" << std::endl;

  Timer t_c;

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
        dash::tasks::async("TRSM",
          [&matrix, k, i, block_size/*, block_b*/](value_t *block_k) {
            Block block_b(matrix, k, i);
            trsm(block_k, block_b.lbegin(), block_size, block_size);
          },
          //DART_PRIO_HIGH,
          (dart_task_prio_t)((num_blocks - i) * (num_blocks-i) * (num_blocks - i) + 3 * ((2 * num_blocks) - k - i - 1) * (i - k))/*priority*/,
          dash::tasks::copyin_r(block_k.begin(), block_k.end()),
          dash::tasks::out(block_b)
        );
        ++num_tasks;
      }
    }
    dash::tasks::async_fence();

    // walk to the right
    for (int i = k+1; i < num_blocks; ++i) {
      Block block_a(matrix, k, i);

      // update diagonal blocks
      Block block_i(matrix, i, i);
      if (block_i.is_local()) {
        // A[j,j] = A[j,j] - A[j,i] * (A[j,i])^t
        dash::tasks::async("SYRK",
          [&matrix, i, block_size/*, block_i*/](value_t *block_a) mutable {
            Block block_i(matrix, i, i);
            syrk(block_a, block_i.lbegin(), block_size, block_size);
          },
          //DART_PRIO_HIGH,
          (dart_task_prio_t)((num_blocks - i) * (num_blocks - i) * (num_blocks - i) + 3 * (i - k))/*priority*/,
          dash::tasks::copyin_r(block_a.begin(), block_a.end()),
          dash::tasks::out(block_i)
        );
        ++num_tasks;
      }

      // run down to the diagonal
      for (int j = i+1; j < num_blocks; ++j) {
        Block block_c(matrix, i, j);
        if (block_c.is_local()) {
          Block block_b(matrix, k, j);
          // A[k,i] = A[k,i] - A[k,j] * (A[j,i])^t
          dash::tasks::async("GEMM",
            [&matrix, block_size, i, j/*, block_c*/](value_t *block_a, value_t *block_b) mutable {
              Block block_c(matrix, i, j);
              gemm(block_a, block_b,
                   block_c.lbegin(), block_size, block_size);
            },
            (dart_task_prio_t)((num_blocks - i) * (num_blocks - i) * (num_blocks - i) + 3 * ((2 * num_blocks) - i - j - 3) * (i - j) + 6 * (i - k)) /*priority*/,
            dash::tasks::copyin_r(block_a.begin(), block_a.end()),
            dash::tasks::copyin_r(block_b.begin(), block_b.end()),
            dash::tasks::out(block_c)
          );
          ++num_tasks;
        }
      }
    }
    dash::tasks::async_fence();
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

}

#endif //CHOLESKY_TASKS_PREFETCH__H
