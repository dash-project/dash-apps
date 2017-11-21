#include <libdash.h>
#include <cblas.h>
#include <lapacke.h>

#include <iostream>
#include <stdlib.h>

using TiledMatrix = dash::Matrix<double, 2, dash::default_index_t, typename dash::ShiftTilePattern<2>>;
static const size_t block_size = 5;


static void fill_random(TiledMatrix &matrix);
static void print_matrix(TiledMatrix &matrix);

class LocalBlockCache {
public:
  LocalBlockCache(
    TiledMatrix &matrix,
    size_t block_row_idx,
    size_t block_col_idx)
  : _matrix(matrix) {
    auto& pattern       = matrix.pattern();
    auto  blockspec     = pattern.blockspec();
    auto  block_extents = blockspec.extents();
    auto  glob_index    = block_row_idx*pattern.blocksize(0)*pattern.extent(1) + block_col_idx*pattern.blocksize(1);
    this->_glob_idx     = glob_index;
    this->_is_local = pattern.is_local(glob_index);
    if (_is_local) {
      _local_ptr = matrix.lbegin() +
                    matrix.pattern().local_index(
                      {block_row_idx*pattern.blocksize(0),
                        block_col_idx*pattern.blocksize(1)}).index;
    }
    this->_size = block_extents[0] * block_extents[1];
    std::cout << "BlockCache: " << block_row_idx << "x" << block_col_idx << " "
              << block_extents[0] << "x" << block_extents[1]
              << " (size=" << _size << ", is_local=" << _is_local
              << ", glob_idx=" << glob_index
              << std::endl;
  }

  ~LocalBlockCache() {
    if (!this->_is_local && _local_ptr != nullptr) {
      std::cout << "Freeing pointer " << _local_ptr << std::endl;
      free(_local_ptr);
      _local_ptr = nullptr;
    }
  }

  double *lbegin() {
    if (!this->_is_local) {
      fetch_data();
    }
    return this->_local_ptr;
  }

  double *lend() {
    if (!this->_is_local) {
      fetch_data();
      return this->_local_ptr + this->_size;
    }
    return this->_local_ptr + this->_size;
  }

  bool is_local() {
    return this->_is_local;
  }

  size_t size() {
    return this->_size;
  }

private:

  void fetch_data() {
    if (this->_local_ptr == nullptr) {
      this->_local_ptr = static_cast<double*>(malloc(this->_size * sizeof(double)));
      auto begin = this->_matrix.begin() + _glob_idx;
      auto end   = begin + _size;
      std::cout << "Fetching " << dash::distance(begin, end) << " elements into "
                << _local_ptr << std::endl;
      dart_get_blocking(
        this->_local_ptr, begin.dart_gptr(),
        dash::distance(begin, end), DART_TYPE_DOUBLE);
      //dash::copy(begin, end, this->_local_ptr);
    }
  }

private:
  TiledMatrix& _matrix;
  double* _local_ptr = nullptr;
  size_t  _size = 0;
  size_t  _glob_idx;
  bool    _is_local  = true;
};

static void
dgemm_tiles(
  TiledMatrix &matrix,
  size_t num_blocks,
  size_t block_size,
  size_t row)
{
  for (size_t j = 0; j < row; ++j) {
    LocalBlockCache block_b(matrix, row, j);
    for (size_t k = row+1; k < num_blocks; ++k) {
      LocalBlockCache block_c(matrix, k, row);
      if (block_c.is_local()) {
        LocalBlockCache block_a(matrix, k, j);
        std::cout << "Calling cblas_dgemm on "
            << block_a.lbegin() << " "
            << block_b.lbegin() << " "
            << block_c.lbegin() << std::endl;
        cblas_dgemm(
            CblasRowMajor, CblasNoTrans, CblasTrans,
            block_size, block_size, block_size,
            -1.0, block_a.lbegin(), block_size,
                  block_b.lbegin(), block_size,
             1.0, block_c.lbegin(), block_size);
      }
    }
  }
}


int main(int argc, char **argv)
{
  dash::init(&argc, &argv);
  const size_t  N = 10 * dash::size();

  TiledMatrix matrix(N, N, dash::TILE(block_size), dash::TILE(block_size));

  fill_random(matrix);

  std::cout << "Matrix: lbegin=" << matrix.lbegin() << ", lend=" << matrix.lend() << std::endl;

  matrix.barrier();
  auto& pattern = matrix.pattern();

  const size_t num_blocks = pattern.blockspec().extent(0);

  if (dash::myid() == 0) {
    std::cout << "block sizes: " << pattern.blocksize(0) << "x" << pattern.blocksize(1) << std::endl;
    std::cout << "num blocks: "  << pattern.blockspec().extent(0) << "x" << pattern.blockspec().extent(1) << std::endl;
  }

  print_matrix(matrix);

  // iterate over column of blocks
  for (size_t i = 0; i < num_blocks; ++i) {

    if (dash::myid() == 0)
      std::cout << "Processing column " << i << std::endl;

    dgemm_tiles(matrix, num_blocks, block_size, i);
    dash::barrier();

    for (size_t j = 0; j < i; ++j) {
      LocalBlockCache block_b(matrix, i, i);
      if (block_b.is_local()) {
        LocalBlockCache block_a(matrix, i, j);
        cblas_dsyrk(
            CblasRowMajor, CblasLower, CblasNoTrans,
            block_size, block_size,
            -1.0, block_a.lbegin(), block_size,
             1.0, block_b.lbegin(), block_size);
      }
    }
    dash::barrier();

    /**
     * Compute cholesky factorization of block on diagonal
     */
//    std::array<typename TiledMatrix::index_type, 2> block_idx = {i, i};
//    auto block = matrix.block(block_idx);
    LocalBlockCache block(matrix, i, i);
    if (block.is_local()) {
      std::cout << "[" << dash::myid() << "] dpotrf on block ("
                << i << ", " << i << ")"
                << " [lbegin: " << block.lbegin() << "]"
                << std::endl;
      LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', block_size, block.lbegin(), block_size);
    }

    dash::barrier();
    /**
     * Solve the triangular equation system in the block
     */
    LocalBlockCache block_a(matrix, i, i);
    for (int j = i+1; j < num_blocks; ++j) {
//      std::array<typename TiledMatrix::index_type, 2> block_b_idx = {i, j};
//      auto& block_b = matrix.block(block_b_idx);
      LocalBlockCache block_b(matrix, i, j);
      if (block_b.is_local()) {
        cblas_dtrsm(
          CblasRowMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit,
          block_size, block_size,
          1.0, block_a.lbegin(), block_size,
               block_b.lbegin(), block_size);
      }
    }

    dash::barrier();
  }

  print_matrix(matrix);


  dash::finalize();

  return 0;
}

static
void fill_random(TiledMatrix &matrix)
{
  constexpr int rand_max = 100;
  dash::generate(
    matrix.begin(),
    matrix.end(),
    [&](){ return (rand())%(rand_max); });
}

static
void print_matrix(TiledMatrix &matrix)
{
  if (matrix.extent(0) > 100) return;

  std::cout << std::setw(5);

  if (dash::myid() == 0) {
    for (size_t i = 0; i < matrix.extent(0); ++i) {
      for (size_t j = 0; j < matrix.extent(1); ++j) {
        std::cout << std::setw(5) << static_cast<double>(matrix(i, j)) << " ";
      }
      std::cout << std::endl;
    }
  }
}

