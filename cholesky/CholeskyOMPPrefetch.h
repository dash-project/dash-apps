#ifndef CHOLESKY_TASKS_PREFETCH__H
#define CHOLESKY_TASKS_PREFETCH__H

#include <libdash.h>
#include <memory>
#include <omp.h>
#include "MatrixBlock.h"
#include "common.h"
#include "ExtraeInstrumentation.h"

constexpr const char *CHOLESKY_IMPL = "CholeskyTasksPrefetch";

template<typename TiledMatrix>
void
compute(TiledMatrix& matrix, size_t block_size){

  using value_t       = typename TiledMatrix::value_type;
  using Block         = MatrixBlock<TiledMatrix>;
  using BlockCache    = typename std::vector<value_t>;
  using BlockCachePtr = typename std::shared_ptr<BlockCache>;
  const size_t num_blocks = matrix.pattern().extent(0) / block_size;

#ifdef USE_EXTRAE
  unsigned nvalues = 6;
  Extrae_define_event_type(&et, "Operations", &nvalues, ev, extrae_names);
#endif

  // pre-allocate a block for pre-fetching the result of potrf()
  value_t *block_k_pre = new value_t[block_size*block_size];
  // allocate a vector to pre-fetch the result of trsm() into
  value_t *blocks_ki_pre = new value_t[block_size*block_size*num_blocks];

  std::cout << "block_k_pre: " << block_k_pre << std::endl;
  std::cout << "blocks_ki_pre: " << blocks_ki_pre << std::endl;
  std::cout << "matrix.lbegin: " << matrix.lbegin() << std::endl;

  /**
   * Algorithm taken from
   * https://pm.bsc.es/ompss-docs/examples/01-examples/cholesky
   */

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
        EXTRAE_ENTER(EVENT_POTRF);
        std::cout << "[" << dash::myid() << ", " << omp_get_thread_num() << "] potrf() on row " << k << "/" << num_blocks << ": ";
        potrf(block_k.lbegin(), block_size, block_size);
        std::cout << "Done." << std::endl;
        EXTRAE_EXIT(EVENT_POTRF);
    }
    dash::barrier();

    // prefetch block_k after it has been computed
    EXTRAE_ENTER(EVENT_PREFETCH);
    dash::dart_storage<value_t> ds(block_size*block_size);
    dart_get_blocking(block_k_pre, block_k.begin().dart_gptr(), ds.nelem, ds.dtype);
    EXTRAE_EXIT(EVENT_PREFETCH);

    std::map<size_t, value_t*> prefetch_blocks;

    /**
     * Solve the triangular equation system in the block
     * Only block_k needs prefetching here.
     */
    for (int i = k+1; i < num_blocks; ++i) {
      Block block_b(matrix, k, i);
      if (block_b.is_local()) {
#pragma omp task firstprivate(block_k_pre, block_b)
{
        EXTRAE_ENTER(EVENT_TRSM);
        trsm(block_k_pre,
            block_b.lbegin(), block_size, block_size);
        EXTRAE_EXIT(EVENT_TRSM);
}
      }
    }
#pragma omp taskwait
    dash::barrier();

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
              auto block_kj_pre = &blocks_ki_pre[j*block_size*block_size];
#pragma omp task depend(out: block_kj_pre[0]) firstprivate(block_kj_pre, block_kj)
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
              auto block_ki_pre = &blocks_ki_pre[i*block_size*block_size];
#pragma omp task depend(out: block_ki_pre[0]) firstprivate(block_ki_pre, block_ki)
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
            auto block_ki_pre = &blocks_ki_pre[i*block_size*block_size];
#pragma omp task depend(out: block_ki_pre[0]) firstprivate(block_ki_pre, block_ki)
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
          auto block_a_pre = prefetch_blocks[i];
          auto block_b_pre = prefetch_blocks[j];
          auto block_c_pre = block_c.lbegin();
          // A[k,i] = A[k,i] - A[k,j] * (A[j,i])^t#
#pragma omp task depend(in: block_a_pre[0], block_b_pre[0]) depend(out: block_c_pre[0]) firstprivate(block_a_pre, block_b_pre, block_c_pre)
{
              EXTRAE_ENTER(EVENT_GEMM);
              gemm(block_a_pre,
                   block_b_pre,
                   block_c_pre, block_size, block_size);
              EXTRAE_EXIT(EVENT_GEMM);
}
        }
      }

      // update diagonal blocks
      Block block_i(matrix, i, i);
      if (block_i.is_local()) {
        //Block block_a(matrix, k, i);
        // A[j,j] = A[j,j] - A[j,i] * (A[j,i])^t
        assert(prefetch_blocks.find(i) != prefetch_blocks.end());
        auto block_a_pre = prefetch_blocks[i];
        auto block_i_pre = block_i.lbegin();
#pragma omp task depend(in: block_a_pre[0], block_i_pre[0]) firstprivate(block_a_pre, block_i_pre)
{
            EXTRAE_ENTER(EVENT_SYRK);
            syrk(block_a_pre,
                 block_i_pre, block_size, block_size);
            EXTRAE_EXIT(EVENT_SYRK);
}
      }
    }
#pragma omp taskwait
    dash::barrier();
  }

  delete[] block_k_pre;
  delete[] blocks_ki_pre;

}

#endif //CHOLESKY_TASKS_PREFETCH__H
