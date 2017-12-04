#ifndef CHOLESKY_TASKS_PREFETCH__H
#define CHOLESKY_TASKS_PREFETCH__H

#include <libdash.h>
#include <memory>
#include "MatrixBlock.h"
#include "common.h"

#ifdef USE_EXTRAE
#include "extrae_user_events.h"
#include "extrae_types.h"
//extern "C" void Extrae_init (void) __attribute__((weak));
//extern "C" void Extrae_event (extrae_type_t type, extrae_value_t value) __attribute__((weak));
//extern "C" void Extrae_fini (void) __attribute__((weak));
#endif 

#define EVENT_POTRF 1000
#define EVENT_TRSM  1001
#define EVENT_GEMM  1002
#define EVENT_SYRK  1003
#define EVENT_PREFETCH 2000

#ifdef USE_EXTRAE
#define EXTRAE_ENTER(_e) Extrae_event(_e, 1)
#define EXTRAE_EXIT(_e)  Extrae_event(_e, 0)
#else
#define EXTRAE_ENTER(_e)
#define EXTRAE_EXIT(_e)
#endif

constexpr const char *CHOLESKY_IMPL = "CholeskyTasksPrefetch";

template<typename TiledMatrix>
void
compute(TiledMatrix& matrix, size_t block_size){

  using value_t       = typename TiledMatrix::value_type;
  using Block         = MatrixBlock<TiledMatrix>;
  using BlockCache    = typename std::vector<value_t>;
  using BlockCachePtr = typename std::shared_ptr<BlockCache>;
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

    auto block_k_pre = std::make_shared<BlockCache>();

    // prefetch block_k after it has been computed
    dash::tasks::async(
      [=]() mutable {
        EXTRAE_ENTER(EVENT_PREFETCH);
        block_k_pre->resize(block_k.size());
        block_k.fetch_async(block_k_pre->data());
        while (!block_k.test()) {
          EXTRAE_EXIT(EVENT_PREFETCH);
          dash::tasks::yield(1);
          EXTRAE_ENTER(EVENT_PREFETCH);
        }
        
        EXTRAE_EXIT(EVENT_PREFETCH);
      },
      //DART_PRIO_HIGH,
      dash::tasks::in(block_k),
      dash::tasks::out(block_k_pre.get())
    );
    ++num_tasks;

    std::map<std::pair<size_t, size_t>, BlockCachePtr> prefetch_blocks;

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
            trsm(block_k_pre->data(),
                block_b.lbegin(), block_size, block_size);
            EXTRAE_EXIT(EVENT_TRSM);

          },
          //DART_PRIO_HIGH,
          dash::tasks::in(block_k_pre.get()),
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
      for (size_t j = k+1; j < i; ++j) {
        Block block_c(matrix, j, i);
        if (block_c.is_local()) {
          auto pkj = std::make_pair(k, j);
          if (prefetch_blocks.find(pkj) == prefetch_blocks.end()) {
            // prefetch block_kj
            auto block_kj_pre = std::make_shared<BlockCache>();
            Block block_kj(matrix, k, j);
            dash::tasks::async(
              [=]() mutable {
                EXTRAE_ENTER(EVENT_PREFETCH);
                block_kj_pre->resize(block_kj.size());
                block_kj.fetch_async(block_kj_pre->data());
                while (!block_kj.test()) {
                  EXTRAE_EXIT(EVENT_PREFETCH);
                  dash::tasks::yield(5);
                  EXTRAE_ENTER(EVENT_PREFETCH);
                }
                EXTRAE_EXIT(EVENT_PREFETCH);
              },
              //DART_PRIO_HIGH,
              dash::tasks::in(block_kj),
              dash::tasks::out(block_kj_pre.get())
            );
            ++num_tasks;
            prefetch_blocks.insert(std::make_pair(pkj, block_kj_pre));
          }
        }
      }
    }


    // walk to the right
    for (size_t i = k+1; i < num_blocks; ++i) {

      auto block_a_pre = std::make_shared<BlockCache>();
      Block block_a(matrix, k, i);
      dash::tasks::async(
        [=]() mutable {
          EXTRAE_ENTER(EVENT_PREFETCH);
          block_a_pre->resize(block_a.size());
          block_a.fetch_async(block_a_pre->data());
          while (!block_a.test()) {
            EXTRAE_EXIT(EVENT_PREFETCH);
            dash::tasks::yield(5);
            EXTRAE_ENTER(EVENT_PREFETCH);
          }
          EXTRAE_EXIT(EVENT_PREFETCH);
        },
        //DART_PRIO_HIGH,
        dash::tasks::in(block_a),
        dash::tasks::out(block_a_pre.get())
      );
      ++num_tasks;

      // run down to the diagonal
      for (size_t j = k+1; j < i; ++j) {
        Block block_c(matrix, j, i);
        if (block_c.is_local()) {
          //Block block_b(matrix, k, j);
          auto block_b_pre = prefetch_blocks[std::make_pair(k, j)];
          // A[k,i] = A[k,i] - A[k,j] * (A[j,i])^t
          dash::tasks::async(
            [=]() mutable {
              EXTRAE_ENTER(EVENT_GEMM);
              gemm(block_a_pre->data(),
                   block_b_pre->data(),
                   block_c.lbegin(), block_size, block_size);
              EXTRAE_EXIT(EVENT_GEMM);
            },
            dash::tasks::in(block_a_pre.get()),
            dash::tasks::in(block_b_pre.get()),
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
        dash::tasks::async(
          [=]() mutable {
            EXTRAE_ENTER(EVENT_SYRK);
            syrk(block_a_pre->data(),
                 block_i.lbegin(), block_size, block_size);
            EXTRAE_EXIT(EVENT_SYRK);
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
