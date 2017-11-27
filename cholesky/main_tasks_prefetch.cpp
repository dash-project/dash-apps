#include <libdash.h>

#include <iostream>
#include <stdlib.h>
#include <dash/util/Timer.h>

#include "common.h"
#include "MatrixBlock.h"

using value_t = double;
using PatternT = typename dash::ShiftTilePattern<2>;
//using PatternT = typename dash::TilePattern<2>;
using TiledMatrix = dash::Matrix<value_t, 2, dash::default_index_t, PatternT>;
using Block = MatrixBlock<TiledMatrix>;
static const size_t N = 2000;
#define BLOCKS_PER_UNIT 10
//#define NUM_BLOCKS      10

//#define DEBUG
//#define CHECK_RESULT

static
void fill_random(TiledMatrix &matrix);
static
void print_matrix(Block &block, size_t nx, size_t ny);

typedef dash::util::Timer<dash::util::TimeMeasure::Clock> Timer;

using BlockCache = typename std::vector<value_t>;

void
compute(TiledMatrix& matrix, size_t block_size){

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

int main(int argc, char **argv)
{
  dash::init(&argc, &argv);
  //const size_t  N = 15 * dash::size();
#ifdef NUM_BLOCKS
  static size_t block_size = N / NUM_BLOCKS;
#elif defined(BLOCKS_PER_UNIT)
  static size_t block_size = N / dash::size() / BLOCKS_PER_UNIT;
#endif


  Timer::Calibrate();

  /*
  if (block_size * dash::size() * BLOCKS_PER_UNIT != N) {
    if (dash::myid() == 0) {
      std::cerr << "Invalid configuration detected: "
                << "  BLOCKS_PER_UNIT: " << BLOCKS_PER_UNIT
                << "  dash::size()   : " << dash::size()
                << "  block_size     : " << block_size
                << std::endl;
    }
    dart_abort(-1);
    return -1;
  }
  */

  TiledMatrix matrix(N, N, dash::TILE(block_size), dash::TILE(block_size));

  fill_random(matrix);

  std::cout << "Matrix: lbegin=" << matrix.lbegin()
            << ", lend=" << matrix.lend() << std::endl;

  matrix.barrier();
  auto& pattern = matrix.pattern();


  if (dash::myid() == 0) {
    std::cout << "block sizes: "
              << pattern.blocksize(0) << "x" << pattern.blocksize(1) << std::endl;
    std::cout << "num blocks: "
              << pattern.blockspec().extent(0) << "x"
              << pattern.blockspec().extent(1) << std::endl;
  }

#ifdef DEBUG
  const size_t num_blocks = pattern.blockspec().extent(0);
  if (N <= 20) {
    print_matrix(matrix);
    if (dash::myid() == 0) {
      std::cout << "\n\nChecking tiles: " << std::endl;
      for (size_t i = 0; i < num_blocks; ++i) {
        for (size_t j = 0; j < num_blocks; ++j) {
          std::cout << "Block " << i << "x" << j << ":\n";
          LocalBlockCache block(matrix, i, j);
          print_matrix(block, block_size, block_size);
          std::cout << std::endl;
        }
      }
    }
  }
#endif

#if defined(CHECK_RESULT)
  TiledMatrix matrix_single(N, N, dash::TILE(block_size), dash::TILE(block_size));
  // copy the matrix before compute
  std::copy(matrix.lbegin(), matrix.lend(), matrix_single.lbegin());
  dash::barrier();
  if (dash::myid() == 0 && N <=20) {
    if (N <=20) {
      print_matrix(matrix_single);
    }
    std::cout << "Computing verification matrix" << std::endl;
  }
  // compute the correct answer on one unit
  compute_single(matrix_single, block_size);
  if (dash::myid() == 0) {
    std::cout << "Done computing verification matrix" << std::endl;
    if (N <=20) {
      std::cout << "########## Expected Result ###############\n";
      print_matrix(matrix_single);
      std::cout << "##########################################\n";
    }
  }
#endif

  compute(matrix, block_size);

  if (N <= 20)
    print_matrix(matrix);

#if defined(CHECK_RESULT)
  verify_matrix(matrix, matrix_single);
#endif
  dash::finalize();

  return 0;
}

static
void fill_random(TiledMatrix &matrix)
{
#if 0 /*defined(DEBUG) || defined(CHECK_RESULT)*/
  // have unit 0 fill the whole matrix
  //int c = 0;
  if (dash::myid() == 0)
  {
    constexpr int rand_max = 100;
    for (auto it = matrix.begin(); it != matrix.end(); ++it) {
      *it = (rand())%(rand_max);
      //*it = c++;
    }
  }
#else
  constexpr int rand_max = 100;
  for (auto it = matrix.lbegin(); it != matrix.lend(); ++it) {
    *it = (rand())%(rand_max);
  }
#endif
  // below causes invalid write!
  /*
  dash::generate(
    matrix.begin(),
    matrix.end(),
    [&](){ return (rand())%(rand_max); });
  */
}

static
void print_matrix(Block &block, size_t nx, size_t ny)
{
  if (dash::myid() == 0) {
    for (size_t i = 0; i < nx; ++i) {
      for (size_t j = 0; j < ny; ++j) {
        std::cout << std::setw(5) << std::setprecision(3)
                  << static_cast<value_t>(*(block.lbegin() + i*ny+j)) << " ";
      }
      std::cout << std::endl;
    }
  }
}


