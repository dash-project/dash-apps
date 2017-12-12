#ifndef CHOLESKY_TASKS__H
#define CHOLESKY_TASKS__H

#include <libdash.h>
#include "MatrixBlock.h"

constexpr const char *CHOLESKY_IMPL = "CholeskyTasks";

template<typename TiledMatrix>
void
compute(TiledMatrix& matrix, size_t block_size){

  using Block = MatrixBlock<TiledMatrix>;
  const size_t num_blocks = matrix.pattern().extent(0) / block_size;
  /**
   * Algorithm taken from
   * https://pm.bsc.es/ompss-docs/examples/01-examples/cholesky
   */

  size_t num_tasks = 0;
  if (dash::myid() == 0)
    std::cout << "Starting to create tasks" << std::endl;

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

    /**
     * Solve the triangular equation system in the block
     */
    for (int i = k+1; i < num_blocks; ++i) {
      Block block_b(matrix, k, i);
      if (block_b.is_local()) {
        dash::tasks::async(
          [=]() mutable {
            trsm(block_k.lbegin(),
                block_b.lbegin(), block_size, block_size);

          },
          dash::tasks::in(block_k),
          dash::tasks::out(block_b)
        );
        ++num_tasks;
      }
    }
    dash::tasks::async_barrier();

    // walk to the right
    for (size_t i = k+1; i < num_blocks; ++i) {
      // run down to the diagonal
      Block block_a(matrix, k, i);
      for (size_t j = k+1; j < i; ++j) {
        Block block_c(matrix, j, i);
        if (block_c.is_local()) {
          Block block_b(matrix, k, j);
          // A[k,i] = A[k,i] - A[k,j] * (A[j,i])^t
          dash::tasks::async(
            [=]() mutable {
              gemm(block_a.lbegin(),
                   block_b.lbegin(),
                   block_c.lbegin(), block_size, block_size);

            },
            dash::tasks::in(block_a),
            dash::tasks::in(block_b),
            dash::tasks::out(block_c)
          );
          ++num_tasks;
        }
      }
    }
    dash::tasks::async_barrier();

    // update diagonal blocks
    for (size_t i = k+1; i < num_blocks; ++i) {
      Block block_i(matrix, i, i);
      if (block_i.is_local()) {
        Block block_a(matrix, k, i);
        // A[j,j] = A[j,j] - A[j,i] * (A[j,i])^t
        dash::tasks::async(
          [=]() mutable {
            syrk(block_a.lbegin(),
                block_i.lbegin(), block_size, block_size);
          },
          dash::tasks::in(block_a),
          dash::tasks::out(block_i)
        );
        ++num_tasks;
      }
    }
    dash::tasks::async_barrier();
  }

  if (dash::myid() == 0)
    std::cout << "Done creating " << num_tasks << " tasks, executing" << std::endl;
  dash::tasks::complete();
  if (dash::myid() == 0)
    std::cout << "Done executing " << num_tasks << " tasks!" << std::endl;


}


#endif // CHOLESKY_TASKS__H
