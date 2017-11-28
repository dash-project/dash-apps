#ifndef CHOLESKY_TASKS_PREFETCH__H
#define CHOLESKY_TASKS_PREFETCH__H

#include <libdash.h>
#include <memory>
#include "MatrixBlock.h"
#include "common.h"

constexpr const char *CHOLESKY_IMPL = "CholeskyTasksPrefetch";

template<typename TiledMatrix>
void
compute(TiledMatrix& matrix, size_t block_size){

  using value_t = typename TiledMatrix::value_type;
  using Block = MatrixBlock<TiledMatrix>;
  using BlockCache = typename std::vector<value_t>;
  const size_t num_blocks = matrix.pattern().blockspec().extent(0);
  /**
   * Algorithm taken from
   * https://pm.bsc.es/ompss-docs/examples/01-examples/cholesky
   */

  size_t num_tasks = 0;
  if (dash::myid() == 0)
    std::cout << "Starting to create tasks" << std::endl;

  Timer t_c;

  // iterate over column of blocks
  for (size_t k = 0; k < num_blocks; ++k) {

    Block block_k(matrix, k, k);

    /**
     * Compute cholesky factorization of block on diagonal
     */
    if (block_k.is_local()) {
      dash::tasks::async(
        [=]() mutable {
          potrf(block_k.lbegin(), block_size, block_size);
        },
        dash::tasks::out(block_k)
      );
      ++num_tasks;
    }
    dash::tasks::async_barrier();

    auto block_k_pre = std::make_shared<BlockCache>();

    // prefetch block_i after it has been computed
    dash::tasks::async(
      [=]() mutable {
        block_k_pre->resize(block_k.size());
        block_k.fetch_async(block_k_pre->data());
        while (!block_k.test())
          dash::tasks::yield(1);
      },
      DART_PRIO_HIGH,
      dash::tasks::in(block_k),
      dash::tasks::out(block_k_pre.get())
    );


    /**
     * Solve the triangular equation system in the block
     */
    for (int i = k+1; i < num_blocks; ++i) {
      Block block_b(matrix, k, i);
      if (block_b.is_local()) {
        dash::tasks::async(
          [=]() mutable {
            trsm(block_k_pre->data(),
                block_b.lbegin(), block_size, block_size);

          },
          DART_PRIO_HIGH,
          dash::tasks::in(block_k_pre.get()),
          dash::tasks::out(block_b)
        );
        ++num_tasks;
      }
    }
    dash::tasks::async_barrier();

    // walk to the right
    for (size_t i = k+1; i < num_blocks; ++i) {
      auto block_a_pre = std::make_shared<BlockCache>();
      Block block_a(matrix, k, i);

      // prefetch block_i after it has been computed
      dash::tasks::async(
        [=]() mutable {
          block_a_pre->resize(block_a.size());
          block_a.fetch_async(block_a_pre->data());
          while (!block_a.test())
            dash::tasks::yield(1);
        },
        DART_PRIO_HIGH,
        dash::tasks::in(block_a),
        dash::tasks::out(block_a_pre.get())
      );

      // run down to the diagonal
      for (size_t j = k+1; j < i; ++j) {
        Block block_c(matrix, j, i);
        if (block_c.is_local()) {
          Block block_b(matrix, k, j);
          // A[k,i] = A[k,i] - A[k,j] * (A[j,i])^t
          dash::tasks::async(
            [=]() mutable {
              gemm(block_a_pre->data(),
                   block_b.lbegin(),
                   block_c.lbegin(), block_size, block_size);

            },
            dash::tasks::in(block_a_pre.get()),
            dash::tasks::in(block_b),
            dash::tasks::out(block_c)
          );
          ++num_tasks;
        }
      }

      // update diagonal blocks
      Block block_i(matrix, i, i);
      if (block_i.is_local()) {
        Block block_a(matrix, k, i);
        // A[j,j] = A[j,j] - A[j,i] * (A[j,i])^t
        dash::tasks::async(
          [=]() mutable {
            syrk(block_a_pre->data(),
                 block_i.lbegin(), block_size, block_size);
          },
          DART_PRIO_HIGH,
          dash::tasks::in(block_a_pre.get()),
          dash::tasks::out(block_i)
        );
        ++num_tasks;
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


}

#endif //CHOLESKY_TASKS_PREFETCH__H
