#ifndef CHOLESKY_TASKS_PREFETCH__H
#define CHOLESKY_TASKS_PREFETCH__H

#include <libdash.h>
#include <memory>
#include "MatrixBlock.h"
#include "common.h"
#include "BlockPrefetcher.h"

constexpr const char *CHOLESKY_IMPL = "CholeskyTasksPrefetch";

template<typename TiledMatrix>
void
compute(TiledMatrix& matrix, size_t block_size){

  using value_t       = typename TiledMatrix::value_type;
  using Block         = MatrixBlock<TiledMatrix>;
  using BlockCache    = typename std::vector<value_t>;
  using BlockCachePtr = typename std::shared_ptr<BlockCache>;
  const size_t num_blocks = matrix.pattern().extent(0) / block_size;

#if 0
  // pre-allocate a block for pre-fetching the result of potrf()
  value_t *block_k_pre = new value_t[block_size*block_size];
  // allocate a vector to pre-fetch the result of trsm() into
  value_t *blocks_ki_pre = new value_t[block_size*block_size*num_blocks];
#endif // 0

  /**
   * Algorithm taken from
   * https://pm.bsc.es/ompss-docs/examples/01-examples/cholesky
   */

  size_t num_tasks = 0;
  if (dash::myid() == 0)
    std::cout << "Starting to create tasks" << std::endl;

  Timer t_c;

  BlockPrefetcher<TiledMatrix> prefetcher(matrix, block_size, num_blocks);

  // iterate over column of blocks
  for (size_t k = 0; k < num_blocks; ++k) {

    Block block_k(matrix, k, k);

    /**
     * Compute cholesky factorization of block on diagonal
     */
    if (block_k.is_local()) {
      dash::tasks::async(
        [=]() mutable {
#ifdef DEBUG
          std::cout << "[" << dash::myid() << ", " << dart_task_thread_num() << "] potrf() on row " << k << "/" << num_blocks << ": ";
#endif
          potrf(block_k.lbegin(), block_size, block_size);
#ifdef DEBUG
          std::cout << "Done." << std::endl;
#endif
        },
        DART_PRIO_HIGH,
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
      auto block_k_pre = prefetcher.get(k, k);
      if (block_b.is_local()) {
        dash::tasks::async(
          [=]() mutable {
            trsm(block_k_pre,
                block_b.lbegin(), block_size, block_size);
          },
          DART_PRIO_HIGH,
          dash::tasks::in(block_k_pre),
          dash::tasks::out(block_b)
        );
        ++num_tasks;
      }
    }
    dash::tasks::async_fence();

    // walk to the right
    for (size_t i = k+1; i < num_blocks; ++i) {

      // run down to the diagonal
      for (size_t j = k+1; j < i; ++j) {
        Block block_c(matrix, j, i);
        if (block_c.is_local()) {
          auto block_a_pre = prefetcher.get(k, i);
          auto block_b_pre = prefetcher.get(k, j);
          // A[k,i] = A[k,i] - A[k,j] * (A[j,i])^t
          dash::tasks::async(
            [=]() mutable {
              gemm(block_a_pre,
                   block_b_pre,
                   block_c.lbegin(), block_size, block_size);
            },
            dash::tasks::in(block_a_pre),
            dash::tasks::in(block_b_pre),
            dash::tasks::out(block_c)
          );
          ++num_tasks;
        }
      }

      // update diagonal blocks
      Block block_i(matrix, i, i);
      if (block_i.is_local()) {
        //Block block_a(matrix, k, i);
        // A[j,j] = A[j,j] - A[j,i] * (A[j,i])^t
        auto block_a_pre = prefetcher.get(k, i);
        dash::tasks::async(
          [=]() mutable {
            syrk(block_a_pre,
                 block_i.lbegin(), block_size, block_size);
          },
          DART_PRIO_HIGH,
          dash::tasks::in(block_a_pre),
          dash::tasks::out(block_i)
        );
        ++num_tasks;
      }
    }
    dash::tasks::async_fence();

    // we need to clear the prefetcher after each iteration
    prefetcher.clear();
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
