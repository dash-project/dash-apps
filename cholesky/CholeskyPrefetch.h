#ifndef CHOLESKY__H
#define CHOLESKY__H

#include <libdash.h>
#include "MatrixBlock.h"
#include "ExtraeInstrumentation.h"

constexpr const char *CHOLESKY_IMPL = "CholeskyNoTasks";

template<typename TiledMatrix>
void
compute(TiledMatrix& matrix, size_t block_size){

  using value_t       = typename TiledMatrix::value_type;
  using Block         = MatrixBlock<TiledMatrix>;
  using BlockCache    = typename std::vector<value_t>;
  using BlockCachePtr = typename std::shared_ptr<BlockCache>;
  const size_t num_blocks = matrix.pattern().extent(0) / block_size;
  /**
   * Algorithm taken from
   * https://pm.bsc.es/ompss-docs/examples/01-examples/cholesky
   */

#ifdef USE_EXTRAE
  unsigned nvalues = 6;
  Extrae_define_event_type(&et, "Operations", &nvalues, ev, extrae_names);
#endif

  // pre-allocate a block for pre-fetching the result of potrf()
  value_t *block_k_pre = new value_t[block_size*block_size];
  // allocate a vector to pre-fetch the result of trsm() into
  value_t *blocks_ki_pre = new value_t[block_size*block_size*num_blocks];


  // iterate over column of blocks
  for (size_t k = 0; k < num_blocks; ++k) {

    if (dash::myid() == 0)
      std::cout << "Processing column " << k << std::endl;

    Block block_k(matrix, k, k);

    /**
     * Compute cholesky factorization of block on diagonal
     */
    if (block_k.is_local()) {
      EXTRAE_ENTER(EVENT_POTRF);
      potrf(block_k.lbegin(), block_size, block_size);
      EXTRAE_EXIT(EVENT_POTRF);
    }
    dash::barrier();

    EXTRAE_ENTER(EVENT_PREFETCH);
    dash::dart_storage<value_t> ds(block_size*block_size);
    dart_get(block_k_pre, block_k.begin().dart_gptr(), ds.nelem, ds.dtype);
    dart_flush(block_k.begin().dart_gptr());
    EXTRAE_EXIT(EVENT_PREFETCH);

    /**
     * Solve the triangular equation system in the block
     */
    for (int i = k+1; i < num_blocks; ++i) {
      Block block_b(matrix, k, i);
      if (block_b.is_local()) {
        EXTRAE_ENTER(EVENT_TRSM);
        trsm(block_k_pre,
             block_b.lbegin(), block_size, block_size);
        EXTRAE_EXIT(EVENT_TRSM);
      }
    }

    dash::barrier();

    std::map<size_t, value_t*> prefetch_blocks;

    /**
     * Prefetch required blocks
     */
    EXTRAE_ENTER(EVENT_PREFETCH);
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
              dash::dart_storage<value_t> ds(block_size*block_size);
              dart_get(block_kj_pre, block_kj.begin().dart_gptr(),
                       ds.nelem, ds.dtype);
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
              dash::dart_storage<value_t> ds(block_size*block_size);
              dart_get(block_ki_pre, block_ki.begin().dart_gptr(),
                       ds.nelem, ds.dtype);
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
            dash::dart_storage<value_t> ds(block_size*block_size);
            dart_get(block_ki_pre, block_ki.begin().dart_gptr(),
                      ds.nelem, ds.dtype);
            prefetch_blocks.insert(
              std::make_pair(i, block_ki_pre));
          }
        }
      }
    }
    dart_flush_all(matrix.begin().dart_gptr());
    EXTRAE_EXIT(EVENT_PREFETCH);

    // walk to the right
    for (size_t i = k+1; i < num_blocks; ++i) {
      // run down to the diagonal
      for (size_t j = k+1; j < i; ++j) {
        Block block_c(matrix, j, i);
        if (block_c.is_local()) {
          auto block_a_pre = prefetch_blocks[i];
          auto block_b_pre = prefetch_blocks[j];
          // A[k,i] = A[k,i] - A[k,j] * (A[j,i])^t
          EXTRAE_ENTER(EVENT_GEMM);
          gemm(block_a_pre,
               block_b_pre,
               block_c.lbegin(), block_size, block_size);
          EXTRAE_EXIT(EVENT_GEMM);
        }
      }

      Block block_ii(matrix, i, i);
      if (block_ii.is_local()) {
        auto block_a_pre = prefetch_blocks[i];
        // A[j,j] = A[j,j] - A[j,i] * (A[j,i])^t
        EXTRAE_ENTER(EVENT_SYRK);
        syrk(block_a_pre,
             block_ii.lbegin(), block_size, block_size);
        EXTRAE_EXIT(EVENT_SYRK);
      }
    }
    dash::barrier();
  }

}

#endif // CHOLESKY__H
