#ifndef CHOLESKY_TASKS_PREFETCH__H
#define CHOLESKY_TASKS_PREFETCH__H

#include <libdash.h>
#include <memory>
#include <mpi.h>
#include "MatrixBlock.h"
#include "common.h"

constexpr const char *CHOLESKY_IMPL = "CholeskyMPIDTasks";

inline static void wait(MPI_Request *comm_req)
{
    int flag = 0;

    MPI_Test(comm_req, &flag, MPI_STATUS_IGNORE);
    while (1) {
        MPI_Test(comm_req, &flag, MPI_STATUS_IGNORE);
        if (flag) break;
        dash::tasks::yield(-1);
    }
//    MPI_Wait(comm_req, MPI_STATUS_IGNORE);
}


inline void reset_send_flags(char *send_flags)
{
    for (int i = 0; i < dash::size(); i++) send_flags[i] = 0;
}

template<typename TiledMatrix>
inline bool block_is_local(TiledMatrix&& matrix, size_t block_size, size_t x, size_t y)
{
  return matrix.pattern.unit_at({x*block_size, y*block_size}) == dash::myid();
}


template<typename TiledMatrix>
void
compute(TiledMatrix& matrix, size_t block_size){

  using value_t       = typename TiledMatrix::value_type;
  using Block         = MatrixBlock<TiledMatrix>;
  using BlockCache    = typename std::vector<value_t>;
  using BlockCachePtr = typename std::shared_ptr<BlockCache>;
  const size_t num_blocks = matrix.pattern().extent(0) / block_size;

  // a bunch of helper variables to ease transition
  size_t np = dash::size();
  size_t nt = num_blocks;
  size_t ts = block_size;
  size_t mype = dash::myid();
  char *send_flags = (char*)std::calloc(sizeof(char), dash::size());

  // pre-allocate a block for pre-fetching the result of potrf()
  value_t *block_k_pre = new value_t[block_size*block_size];
  // allocate a vector to pre-fetch the result of trsm() into
  value_t *blocks_ki_pre = new value_t[block_size*block_size*num_blocks];

  char recv_flag   = 0;

  double ***A = new double**[nt];
  int *block_rank = new int[nt*nt];
  double *B  = block_k_pre;
  double **C = new double*[nt];
  for (int i = 0; i < nt; ++i) {
    A[i] = new double*[nt];
    for (int j = 0; j < nt; ++j) {
      Block block(matrix, i, j);
      A[i][j] = (block.is_local()) ? block.lbegin() : NULL;
      block_rank[i*nt+j] = block.unit();
    }
    C[i] = &blocks_ki_pre[i*block_size*block_size];
  }

  /**
   * Algorithm taken from
   * https://pm.bsc.es/ompss-docs/examples/
   */

  if (dash::myid() == 0)
    std::cout << "Starting to create tasks" << std::endl;

  // used to make sure that all trsm tasks finished before the MPI task runs
  int trsm_sentinel;
  (void)trsm_sentinel;

  // used to ensure that MPI send/recv task finishes before dgemm/syrk runs
  int mpi_sentinel;
  (void)mpi_sentinel;

  Timer t_c;

    for (int k = 0; k < nt; k++) {
        if (block_rank[k*nt+k] == mype) {
          dash::tasks::async(
            [=](){potrf(A[k][k], ts, ts);},
            dash::tasks::out(A[k][k])
          );
        }

        if (block_rank[k*nt+k] == mype && np != 1) {
            for (int kk = k+1; kk < nt; kk++) {
                if (!send_flags[block_rank[k*nt+kk]]) send_flags[block_rank[k*nt+kk]] = 1;
            }
            for (int dst = 0; dst < np; dst++) {
                if (send_flags[dst] && dst != mype) {
                  dash::tasks::async(
                    [=, &A](){
                      MPI_Request send_req;
                      MPI_Isend(A[k][k], ts*ts, MPI_DOUBLE, dst, k*nt+k, MPI_COMM_WORLD, &send_req);
                      wait(&send_req);
                    },
                    dash::tasks::in(A[k][k])
                  );
                }
            }
            reset_send_flags(send_flags);
        }

        if (block_rank[k*nt+k] != mype) {
            for (int i = k + 1; i < nt; i++) {
                if (block_rank[k*nt+i] == mype) recv_flag = 1;
            }
            if (recv_flag) {
              dash::tasks::async(
                [=](){
                  MPI_Request recv_req;
                  MPI_Irecv(B, ts*ts, MPI_DOUBLE, block_rank[k*nt+k], k*nt+k, MPI_COMM_WORLD, &recv_req);
                  wait(&recv_req);
                },
                dash::tasks::out(B)
              );
                recv_flag = 0;
            }
        }

        for (int i = k + 1; i < nt; i++) {
            if (block_rank[k*nt+i] == mype) {
                if (block_rank[k*nt+k] == mype) {
                  dash::tasks::async(
                    [=](){trsm(A[k][k], A[k][i], ts, ts);},
                    dash::tasks::in(A[k][k]),
                    dash::tasks::out(A[k][i])
                  );
                } else {
                  dash::tasks::async(
                    [=](){trsm(B, A[k][i], ts, ts);},
                    dash::tasks::in(B),
                    dash::tasks::out(A[k][i])
                  );
                }
            }

            if (block_rank[k*nt+i] == mype && np != 1) {
                for (int ii = k + 1; ii < i; ii++) {
                    if (!send_flags[block_rank[ii*nt+i]]) send_flags[block_rank[ii*nt+i]] = 1;
                }
                for (int ii = i + 1; ii < nt; ii++) {
                    if (!send_flags[block_rank[i*nt+ii]]) send_flags[block_rank[i*nt+ii]] = 1;
                }
                if (!send_flags[block_rank[i*nt+i]]) send_flags[block_rank[i*nt+i]] = 1;
                for (int dst = 0; dst < np; dst++) {
                    if (send_flags[dst] && dst != mype) {
                      dash::tasks::async(
                        [=](){
                          MPI_Request send_req;
                          MPI_Isend(A[k][i], ts*ts, MPI_DOUBLE, dst, k*nt+i, MPI_COMM_WORLD, &send_req);
                          wait(&send_req);
                        },
                        dash::tasks::in(A[k][i])
                      );
                    }
                }
                reset_send_flags(send_flags);
            }
            if (block_rank[k*nt+i] != mype) {
                for (int ii = k + 1; ii < i; ii++) {
                    if (block_rank[ii*nt+i] == mype) recv_flag = 1;
                }
                for (int ii = i + 1; ii < nt; ii++) {
                    if (block_rank[i*nt+ii] == mype) recv_flag = 1;
                }
                if (block_rank[i*nt+i] == mype) recv_flag = 1;
                if (recv_flag) {
                  dash::tasks::async(
                    [=](){
                      MPI_Request recv_req;
                      MPI_Irecv(C[i], ts*ts, MPI_DOUBLE, block_rank[k*nt+i], k*nt+i, MPI_COMM_WORLD, &recv_req);
                      wait(&recv_req);
                    },
                    dash::tasks::out(C[i])
                  );
                  recv_flag = 0;
                }
            }
        }

        for (int i = k + 1; i < nt; i++) {

            for (int j = k + 1; j < i; j++) {
                if (block_rank[j*nt+i] == mype) {
                    if (block_rank[k*nt+i] == mype && block_rank[k*nt+j] == mype) {
                      dash::tasks::async(
                        [=](){ gemm(A[k][i], A[k][j], A[j][i], ts, ts); },
                        dash::tasks::in(A[k][i]),
                        dash::tasks::in(A[k][j]),
                        dash::tasks::out(A[j][i])
                      );
                    } else if (block_rank[k*nt+i] != mype && block_rank[k*nt+j] == mype) {
                      dash::tasks::async(
                        [=](){ gemm(C[i], A[k][j], A[j][i], ts, ts); },
                        dash::tasks::in(C[i]),
                        dash::tasks::in(A[k][j]),
                        dash::tasks::out(A[j][i])
                      );
                    } else if (block_rank[k*nt+i] == mype && block_rank[k*nt+j] != mype) {
                      dash::tasks::async(
                        [=](){ gemm(A[k][i], C[j], A[j][i], ts, ts); },
                        dash::tasks::in(A[k][i]),
                        dash::tasks::in(C[j]),
                        dash::tasks::out(A[j][i])
                      );
                    } else {
                      dash::tasks::async(
                        [=](){ gemm(C[i], C[j], A[j][i], ts, ts); },
                        dash::tasks::in(C[i]),
                        dash::tasks::in(C[j]),
                        dash::tasks::out(A[j][i])
                      );
                    }
                }
            }

            if (block_rank[i*nt+i] == mype) {
                if (block_rank[k*nt+i] == mype) {
                  dash::tasks::async(
                    [=](){ syrk(A[k][i], A[i][i], ts, ts); },
                    dash::tasks::in(A[k][i]),
                    dash::tasks::out(A[i][i])
                  );
                } else {
                  dash::tasks::async(
                    [=](){ syrk(C[i], A[i][i], ts, ts); },
                    dash::tasks::in(C[i]),
                    dash::tasks::out(A[i][i])
                  );
                }
            }
        }
    }
  if (dash::myid() == 0)
    std::cout << "Done creating tasks in "
              << t_c.Elapsed() / 1E3 << "ms"
              << ", starting execution" << std::endl;
  dash::tasks::complete();
  if (dash::myid() == 0)
    std::cout << "Done executing tasks in "
              << t_c.Elapsed() / 1E3 << "ms"
              << std::endl;
  delete[] block_k_pre;
  delete[] blocks_ki_pre;

  for (int i = 0; i < nt; ++i) {
    delete[] A[i];
  }
  delete[] A;
  delete[] C;
  delete[] block_rank;
  std::free(send_flags);

}

#endif //CHOLESKY_TASKS_PREFETCH__H
