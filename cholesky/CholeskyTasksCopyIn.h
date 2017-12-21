#ifndef CHOLESKY_TASKS_PREFETCH__H
#define CHOLESKY_TASKS_PREFETCH__H

#include <libdash.h>
#include <memory>
#include "MatrixBlock.h"
#include "common.h"
#include "ExtraeInstrumentation.h"
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

  std::cout << "local_size: " << matrix.local_size() << std::endl;


  // allocate a vector to pre-fetch the result of trsm() into
  value_t *blocks_ki_pre = new value_t[block_size*block_size*num_blocks];

  auto next_buf_pos = [=](){
    static size_t buffer_pos = 0;
    size_t res = buffer_pos;
    buffer_pos = (buffer_pos + 1) % num_blocks;
    return res*block_size*block_size;
  };


#ifdef USE_EXTRAE
  unsigned nvalues = 6;
  Extrae_define_event_type(&et, "Operations", &nvalues, ev, extrae_names);
#endif

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

  // iterate over column of blocks
  for (size_t k = 0; k < num_blocks; ++k) {

    Block block_k(matrix, k, k);

    /**
     * Compute cholesky factorization of block on diagonal
     */
    if (block_k.is_local()) {
      dash::tasks::async(
        [=]() mutable {
          EXTRAE_ENTER(EVENT_POTRF);
          std::cout << "[" << dash::myid() << ", " << dart_task_thread_num()
                    << "] potrf() on row " << k << "/" << num_blocks << ": ";
          potrf(block_k.lbegin(), block_size, block_size);
          std::cout << "Done." << std::endl;
          EXTRAE_EXIT(EVENT_POTRF);
        },
        DART_PRIO_HIGH,
        dash::tasks::out(block_k)
      );
      ++num_tasks;
    }
    dash::tasks::async_barrier();

    /**
     * Solve the triangular equation system in the block
     * Only block_k needs prefetching here.
     */
    for (int i = k+1; i < num_blocks; ++i) {
      Block block_b(matrix, k, i);
      if (block_b.is_local()) {
        //auto block_k_pre = &blocks_ki_pre[k*block_size*block_size];
        auto block_k_pre = block_k.is_local()
                              ? block_k.lbegin()
                              : &blocks_ki_pre[k*block_size*block_size];
        auto block_a_dep = block_k.is_local()
                              ? dash::tasks::in(block_k)
                              : dash::tasks::copyin(block_k, block_size*block_size, block_k_pre);
        dash::tasks::async(
          [=]() mutable {
            EXTRAE_ENTER(EVENT_TRSM);
            trsm(block_k_pre,
                block_b.lbegin(), block_size, block_size);
            EXTRAE_EXIT(EVENT_TRSM);
          },
          DART_PRIO_HIGH,
          block_a_dep,
          dash::tasks::out(block_b)
        );
        ++num_tasks;
      }
    }
    dash::tasks::async_barrier();

    // walk to the right
    for (size_t i = k+1; i < num_blocks; ++i) {
      Block block_a(matrix, k, i);
      auto block_a_pre = block_a.is_local()
                            ? block_a.lbegin()
                            : &blocks_ki_pre[i*block_size*block_size];
      auto block_a_dep = block_a.is_local()
                            ? dash::tasks::in(block_a)
                            : dash::tasks::copyin(block_a, block_size*block_size, block_a_pre);
      // run down to the diagonal
      for (size_t j = k+1; j < i; ++j) {
        Block block_c(matrix, j, i);
        if (block_c.is_local()) {
          Block block_b(matrix, k, j);
          auto block_b_pre = block_b.is_local()
                              ? block_b.lbegin()
                              : &blocks_ki_pre[j*block_size*block_size];
          auto block_b_dep = block_b.is_local()
                              ? dash::tasks::in(block_b)
                              : dash::tasks::copyin(block_b, block_size*block_size, block_b_pre);
          // A[k,i] = A[k,i] - A[k,j] * (A[j,i])^t
          dash::tasks::async(
            [=]() mutable {
              EXTRAE_ENTER(EVENT_GEMM);
              gemm(block_a_pre,
                   block_b_pre,
                   block_c.lbegin(), block_size, block_size);
              EXTRAE_EXIT(EVENT_GEMM);
            },
            block_a_dep,
            block_b_dep,
            dash::tasks::out(block_c)
          );
          ++num_tasks;
        }
      }

      // update diagonal blocks
      Block block_i(matrix, i, i);
      if (block_i.is_local()) {
        // A[j,j] = A[j,j] - A[j,i] * (A[j,i])^t
        dash::tasks::async(
          [=]() mutable {
            EXTRAE_ENTER(EVENT_SYRK);
            syrk(block_a_pre,
                 block_i.lbegin(), block_size, block_size);
            EXTRAE_EXIT(EVENT_SYRK);
          },
          DART_PRIO_HIGH,
          block_a_dep,
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

  delete[] blocks_ki_pre;

}

#endif //CHOLESKY_TASKS_PREFETCH__H
