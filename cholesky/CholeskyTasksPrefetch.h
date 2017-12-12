#ifndef CHOLESKY_TASKS_PREFETCH__H
#define CHOLESKY_TASKS_PREFETCH__H

#include <libdash.h>
#include <memory>
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

  std::cout << "local_size: " << matrix.local_size() << std::endl;

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
          std::cout << "[" << dash::myid() << ", " << dart_task_thread_num() << "] potrf() on row " << k << "/" << num_blocks << ": ";
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

    // prefetch block_k after it has been computed
    dash::tasks::async(
      [=]() mutable {
        EXTRAE_ENTER(EVENT_PREFETCH);
        block_k.fetch_async(block_k_pre);
        while (!block_k.test()) {
          EXTRAE_EXIT(EVENT_PREFETCH);
          dash::tasks::yield(-1);
          EXTRAE_ENTER(EVENT_PREFETCH);
        }

        EXTRAE_EXIT(EVENT_PREFETCH);
      },
      DART_PRIO_HIGH,
      dash::tasks::in(block_k),
      dash::tasks::out(block_k_pre)
    );
    ++num_tasks;

    std::map<size_t, value_t*> prefetch_blocks;

    /**
     * Solve the triangular equation system in the block
     * Only block_k needs prefetching here.
     */
    for (int i = k+1; i < num_blocks; ++i) {
      Block block_b(matrix, k, i);
      if (block_b.is_local()) {
        dash::tasks::async(
          [=]() mutable {
            EXTRAE_ENTER(EVENT_TRSM);
            trsm(block_k_pre,
                block_b.lbegin(), block_size, block_size);
            EXTRAE_EXIT(EVENT_TRSM);

          },
          DART_PRIO_HIGH,
          dash::tasks::in(block_k_pre),
          dash::tasks::out(block_b)
        );
        ++num_tasks;
      }
    }
    dash::tasks::async_barrier();

    /**
     * Prefetch required blocks
     */
    for (size_t i = k+1; i < num_blocks; ++i) {

      Block block_ki(matrix, k, i); // result of trsm() above
      Block block_ii(matrix, i, i); // diagonal blocks

      for (size_t j = k+1; j < i; ++j) {
        Block block_c(matrix, j, i); // lower triangular blocks below current row
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
              dash::tasks::async(
                [=]() mutable {
                  EXTRAE_ENTER(EVENT_PREFETCH);
                  block_kj.fetch_async(block_kj_pre);
                  while (!block_kj.test()) {
                    EXTRAE_EXIT(EVENT_PREFETCH);
                    dash::tasks::yield(-1);
                    EXTRAE_ENTER(EVENT_PREFETCH);
                  }
                  EXTRAE_EXIT(EVENT_PREFETCH);
                },
                //DART_PRIO_HIGH,
                dash::tasks::in(block_kj),
                dash::tasks::out(block_kj_pre)
              );
              ++num_tasks;
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
              dash::tasks::async(
                [=]() mutable {
                  EXTRAE_ENTER(EVENT_PREFETCH);
                  block_ki.fetch_async(block_ki_pre);
                  while (!block_ki.test()) {
                    EXTRAE_EXIT(EVENT_PREFETCH);
                    dash::tasks::yield(-1);
                    EXTRAE_ENTER(EVENT_PREFETCH);
                  }
                  EXTRAE_EXIT(EVENT_PREFETCH);
                },
                //DART_PRIO_HIGH,
                dash::tasks::in(block_ki),
                dash::tasks::out(block_ki_pre)
              );
              ++num_tasks;
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
            dash::tasks::async(
              [=]() mutable {
                EXTRAE_ENTER(EVENT_PREFETCH);
                block_ki.fetch_async(block_ki_pre);
                while (!block_ki.test()) {
                  EXTRAE_EXIT(EVENT_PREFETCH);
                  dash::tasks::yield(-1);
                  EXTRAE_ENTER(EVENT_PREFETCH);
                }
                EXTRAE_EXIT(EVENT_PREFETCH);
              },
              //DART_PRIO_HIGH,
              dash::tasks::in(block_ki),
              dash::tasks::out(block_ki_pre)
            );
            ++num_tasks;
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
          // A[k,i] = A[k,i] - A[k,j] * (A[j,i])^t
          dash::tasks::async(
            [=]() mutable {
              EXTRAE_ENTER(EVENT_GEMM);
              gemm(block_a_pre,
                   block_b_pre,
                   block_c.lbegin(), block_size, block_size);
              EXTRAE_EXIT(EVENT_GEMM);
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
        assert(prefetch_blocks.find(i) != prefetch_blocks.end());
        auto block_a_pre = prefetch_blocks[i];
        dash::tasks::async(
          [=]() mutable {
            EXTRAE_ENTER(EVENT_SYRK);
            syrk(block_a_pre,
                 block_i.lbegin(), block_size, block_size);
            EXTRAE_EXIT(EVENT_SYRK);
          },
          DART_PRIO_HIGH,
          dash::tasks::in(block_a_pre),
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

  delete[] block_k_pre;
  delete[] blocks_ki_pre;

}

#endif //CHOLESKY_TASKS_PREFETCH__H
