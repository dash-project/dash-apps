#ifndef CHOLESKY_TASKS_PREFETCH__H
#define CHOLESKY_TASKS_PREFETCH__H

#include <libdash.h>
#include <memory>
#include <omp.h>
#include "MatrixBlock.h"
#include "common.h"

constexpr const char *CHOLESKY_IMPL = "CholeskyOpenMPPrefetch";

template<typename TiledMatrix>
void
compute(TiledMatrix& matrix, size_t block_size){

  using value_t       = typename TiledMatrix::value_type;
  using Block         = MatrixBlock<TiledMatrix>;
  using BlockCache    = typename std::vector<value_t>;
  using BlockCachePtr = typename std::shared_ptr<BlockCache>;
  const size_t num_blocks = matrix.pattern().extent(0) / block_size;

  // pre-allocate a block for pre-fetching the result of potrf()
  value_t *block_k_pre = new value_t[block_size*block_size];
  // allocate a vector to pre-fetch the result of trsm() into
  value_t *blocks_ki_pre = new value_t[block_size*block_size*num_blocks];

  int barrier_sentinel; // dummy variable for barrier synchronization

  /**
   * Algorithm taken from
   * https://pm.bsc.es/ompss-docs/examples/01-examples/cholesky
   */

  if (dash::myid() == 0)
    std::cout << "Starting to create tasks" << std::endl;

  Timer t_c;

#pragma omp parallel shared(matrix, barrier_sentinel)
{
#pragma omp master
{
  // iterate over column of blocks
  for (size_t k = 0; k < num_blocks; ++k) {

    Block block_k(matrix, k, k);

    /**
     * Compute cholesky factorization of block on diagonal
     */
    if (block_k.is_local()) {
      auto block_k_ptr = block_k.lbegin();
#pragma omp task depend(out:block_k_ptr, barrier_sentinel) firstprivate(block_k)
      {
#ifdef DEBUG
        std::cout << "[" << dash::myid() << ", " << omp_get_thread_num() << "] potrf() on row " << k << "/" << num_blocks << ": ";
#endif
        potrf(block_k_ptr, block_size, block_size);
#ifdef DEBUG
        std::cout << "Done." << std::endl;
#endif
      }
    }

#pragma omp task depend(out: barrier_sentinel)
    {dash::barrier();}

    // prefetch block_k after it has been computed
#pragma omp task depend(out: block_k_pre) depend(in: barrier_sentinel) firstprivate(block_k_pre, block_k)
  {
    dash::dart_storage<value_t> ds(block_size*block_size);
    dart_get_blocking(block_k_pre, block_k.begin().dart_gptr(), ds.nelem, ds.dtype, ds.dtype);
  }
    std::map<size_t, value_t*> prefetch_blocks;

    /**
     * Solve the triangular equation system in the block
     * Only block_k needs prefetching here.
     */
    for (int i = k+1; i < num_blocks; ++i) {
      Block block_b(matrix, k, i);
      if (block_b.is_local()) {
        auto block_b_ptr = block_b.lbegin();
#pragma omp task depend(in: block_k_pre, block_b_ptr, barrier_sentinel) firstprivate(block_k_pre, block_b_ptr)
  {
        trsm(block_k_pre,
             block_b_ptr, block_size, block_size);
  }
      }
    }
#pragma omp task depend(out: barrier_sentinel)
    {dash::barrier();}

    /**
     * Prefetch required blocks
     */
    for (size_t i = k+1; i < num_blocks; ++i) {

      Block block_ki(matrix, k, i);
      Block block_ii(matrix, i, i);

      for (size_t j = k+1; j < i; ++j) {
        Block block_c(matrix, j, i);
        if (block_c.is_local()) {
          if (prefetch_blocks.find(j) == prefetch_blocks.end()) {
            // prefetch block_kj
            Block block_kj(matrix, k, j);
            if (block_kj.is_local()) {
              // use local pointer
              prefetch_blocks.insert(
                std::make_pair(j, block_kj.lbegin()));
            } else {
              value_t* block_kj_pre = &blocks_ki_pre[j*block_size*block_size];
#pragma omp task firstprivate(block_kj_pre, block_kj) depend(out: block_kj_pre[0]) depend(in: barrier_sentinel)
{
                  EXTRAE_ENTER(EVENT_PREFETCH);
                  block_kj.fetch_async(block_kj_pre);
                  while (!block_kj.test()) {
                    EXTRAE_EXIT(EVENT_PREFETCH);
                    #pragma omp taskyield
                    EXTRAE_ENTER(EVENT_PREFETCH);
                  }
                  EXTRAE_EXIT(EVENT_PREFETCH);
}
              prefetch_blocks.insert(std::make_pair(j, block_kj_pre));
            }
          }
          if (prefetch_blocks.find(i) == prefetch_blocks.end()) {
            // pre-fetch block_ki
            if (block_ki.is_local()) {
              // use local pointer
              prefetch_blocks.insert(
                std::make_pair(i, block_ki.lbegin()));
            } else {
              // pre-fetch
              value_t* block_ki_pre = &blocks_ki_pre[i*block_size*block_size];
#pragma omp task depend(out: block_ki_pre[0]) firstprivate(block_ki_pre, block_ki) depend(in: barrier_sentinel)
{
                  EXTRAE_ENTER(EVENT_PREFETCH);
                  block_ki.fetch_async(block_ki_pre);
                  while (!block_ki.test()) {
                    EXTRAE_EXIT(EVENT_PREFETCH);
                    #pragma omp taskyield
                    EXTRAE_ENTER(EVENT_PREFETCH);
                  }
                  EXTRAE_EXIT(EVENT_PREFETCH);
}
              prefetch_blocks.insert(
                std::make_pair(i, block_ki_pre));
            }
          }
        }
      }

      // pre-fetch block_ki
      if (block_ii.is_local()) {
        if (prefetch_blocks.find(i) == prefetch_blocks.end()) {
          // pre-fetch block_ki
          if (block_ki.is_local()) {
            // use local pointer
            prefetch_blocks.insert(
              std::make_pair(i, block_ki.lbegin()));
          } else {
            // pre-fetch
            value_t* block_ki_pre = &blocks_ki_pre[i*block_size*block_size];
#pragma omp task depend(out: block_ki_pre[0]) firstprivate(block_ki_pre, block_ki) depend(in: barrier_sentinel)
{
                EXTRAE_ENTER(EVENT_PREFETCH);
                block_ki.fetch_async(block_ki_pre);
                while (!block_ki.test()) {
                  EXTRAE_EXIT(EVENT_PREFETCH);
                    #pragma omp taskyield
                  EXTRAE_ENTER(EVENT_PREFETCH);
                }
                EXTRAE_EXIT(EVENT_PREFETCH);
}
            prefetch_blocks.insert(
              std::make_pair(i, block_ki_pre));
          }
        }
      }
    }


    // walk to the right
    for (size_t i = k+1; i < num_blocks; ++i) {

      // run down to the diagonal
      for (size_t j = k+1; j < i; ++j) {
        Block block_c(matrix, j, i);
        if (block_c.is_local()) {
          //Block block_b(matrix, k, j);
          assert(prefetch_blocks.find(i) != prefetch_blocks.end());
          value_t* block_a_pre = prefetch_blocks[i];
          value_t* block_b_pre = prefetch_blocks[j];
          value_t* block_c_pre = block_c.lbegin();
          // A[k,i] = A[k,i] - A[k,j] * (A[j,i])^t#
#pragma omp task depend(in: block_a_pre[0], block_b_pre[0], barrier_sentinel) depend(out: block_c_pre[0]) firstprivate(block_a_pre, block_b_pre, block_c_pre)
{
              gemm(block_a_pre,
                   block_b_pre,
                   block_c_pre, block_size, block_size);
}
        }
      }

      // update diagonal blocks
      Block block_i(matrix, i, i);
      if (block_i.is_local()) {
        //Block block_a(matrix, k, i);
        // A[j,j] = A[j,j] - A[j,i] * (A[j,i])^t
        assert(prefetch_blocks.find(i) != prefetch_blocks.end());
        value_t* block_a_pre = prefetch_blocks[i];
        value_t* block_i_pre = block_i.lbegin();
#pragma omp task depend(in: block_a_pre[0], block_i_pre[0], barrier_sentinel) firstprivate(block_a_pre, block_i_pre)
{
            syrk(block_a_pre,
                 block_i_pre, block_size, block_size);
}
      }
    }
#pragma omp task depend(out: barrier_sentinel)
    {dash::barrier();}
  }
} // pragma omp master
} // pragma omp parallel

  delete[] block_k_pre;
  delete[] blocks_ki_pre;

}

#endif //CHOLESKY_TASKS_PREFETCH__H
