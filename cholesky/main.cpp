#include <libdash.h>
#include <cblas.h>
#include <lapacke.h>

#include <iostream>
#include <stdlib.h>

#include "common.h"
#include "MatrixBlock.h"

using value_t = double;
using PatternT = typename dash::ShiftTilePattern<2>;
//using PatternT = typename dash::TilePattern<2>;
using TiledMatrix = dash::Matrix<value_t, 2, dash::default_index_t, PatternT>;
using Block = MatrixBlock<TiledMatrix>;
static const size_t N = 20;
#define BLOCKS_PER_UNIT 5
#define NUM_BLOCKS      10

//#define DEBUG
//#define CHECK_RESULT

static void fill_random(TiledMatrix &matrix);
static
void print_matrix(Block &block, size_t nx, size_t ny);

void
compute(TiledMatrix& matrix, size_t block_size){

  const size_t num_blocks = matrix.pattern().blockspec().extent(0);
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

int main(int argc, char **argv)
{
  dash::init(&argc, &argv);
  //const size_t  N = 15 * dash::size();
#ifdef NUM_BLOCKS
  static size_t block_size = N / NUM_BLOCKS;
#elif defined(BLOCKS_PER_UNIT)
  static size_t block_size = N / dash::size() / BLOCKS_PER_UNIT;
#endif

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


  dash::finalize();

  return 0;
}

static
void fill_random(TiledMatrix &matrix)
{
#if defined(DEBUG) || defined(CHECK_RESULT)
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


