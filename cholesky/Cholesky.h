#ifndef CHOLESKY__H
#define CHOLESKY__H

#include <libdash.h>
#include "MatrixBlock.h"

constexpr const char *CHOLESKY_IMPL = "CholeskyNoTasks";

template<typename TiledMatrix>
void
compute(TiledMatrix& matrix, size_t block_size){

  using Block = MatrixBlock<TiledMatrix>;
  const size_t num_blocks = matrix.pattern().extent(0) / block_size;
  /**
   * Algorithm taken from
   * https://pm.bsc.es/ompss-docs/examples/01-examples/cholesky
   */

  // iterate over column of blocks
  for (size_t k = 0; k < num_blocks; ++k) {

    if (dash::myid() == 0)
      std::cout << "Processing column " << k << std::endl;

    Block block_k(matrix, k, k);

    /**
     * Compute cholesky factorization of block on diagonal
     */
    if (block_k.is_local()) {
      potrf(block_k.lbegin(), block_size, block_size);
    }
    dash::barrier();

    /**
     * Solve the triangular equation system in the block
     */
    for (int i = k+1; i < num_blocks; ++i) {
      Block block_b(matrix, k, i);
      if (block_b.is_local()) {
        trsm(block_k.lbegin(),
             block_b.lbegin(), block_size, block_size);
      }
    }

    dash::barrier();

    // walk to the right
    for (size_t i = k+1; i < num_blocks; ++i) {
      // run down to the diagonal
      Block block_a(matrix, k, i);
      for (size_t j = k+1; j < i; ++j) {
        Block block_c(matrix, j, i);
        if (block_c.is_local()) {
          Block block_b(matrix, k, j);
          // A[k,i] = A[k,i] - A[k,j] * (A[j,i])^t
          gemm(block_a.lbegin(),
               block_b.lbegin(),
               block_c.lbegin(), block_size, block_size);
        }
      }
    }
    dash::barrier();

    // update diagonal blocks
    for (size_t i = k+1; i < num_blocks; ++i) {
      Block block_i(matrix, i, i);
      if (block_i.is_local()) {
        Block block_a(matrix, k, i);
        // A[j,j] = A[j,j] - A[j,i] * (A[j,i])^t
        syrk(block_a.lbegin(),
             block_i.lbegin(), block_size, block_size);
      }
    }
    dash::barrier();
  }

}

#endif // CHOLESKY__H
