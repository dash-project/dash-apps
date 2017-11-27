#include <libdash.h>
#include <cblas.h>
#include <lapacke.h>

#include <iostream>
#include <stdlib.h>

using value_t = double;
using PatternT = typename dash::ShiftTilePattern<2>;
//using PatternT = typename dash::TilePattern<2>;
using TiledMatrix = dash::Matrix<value_t, 2, dash::default_index_t, PatternT>;
static const size_t N = 20;
#define BLOCKS_PER_UNIT 5
#define NUM_BLOCKS      10

//#define DEBUG
#define CHECK_RESULT

extern "C" {
void dgemm_ (const char *transa, const char *transb, int *l, int *n, int *m, double *alpha,
             const void *a, int *lda, void *b, int *ldb, double *beta, void *c, int *ldc);
void dtrsm_ (char *side, char *uplo, char *transa, char *diag, int *m, int *n, double *alpha,
             double *a, int *lda, double *b, int *ldb);
void dtrmm_ (char *side, char *uplo, char *transa, char *diag, int *m, int *n, double *alpha,
             double *a, int *lda, double *b, int *ldb);
void dsyrk_ (char *uplo, char *trans, int *n, int *k, double *alpha, double *a, int *lda,
             double *beta, double *c, int *ldc);
//void dpotrf_(const char* uplo, int& n, float* a, int* lda, int* info);
}

void potrf(double * const A, int ts, int ld)
{
   static int INFO;
   static char L = 'L';
   dpotrf_(&L, &ts, A, &ld, &INFO);
}

void trsm(double *A, double *B, int ts, int ld)
{
   static char LO = 'L', TR = 'T', NU = 'N', RI = 'R';
   static double DONE = 1.0;
   dtrsm_(&RI, &LO, &TR, &NU, &ts, &ts, &DONE, A, &ld, B, &ld );
}

void syrk(double *A, double *B, int ts, int ld)
{
   static char LO = 'L', NT = 'N';
   static double DONE = 1.0, DMONE = -1.0;
   dsyrk_(&LO, &NT, &ts, &ts, &DMONE, A, &ld, &DONE, B, &ld );
}

void gemm(double *A, double *B, double *C, int ts, int ld)
{
   static const char TR = 'T', NT = 'N';
   static double DONE = 1.0, DMONE = -1.0;
   dgemm_(&NT, &TR, &ts, &ts, &ts, &DMONE, A, &ld, B, &ld, &DONE, C, &ld);
}

static void fill_random(TiledMatrix &matrix);
static void print_matrix(TiledMatrix &matrix);
class LocalBlockCache;
static
void print_matrix(LocalBlockCache &block, size_t nx, size_t ny);

class LocalBlockCache {
public:
  LocalBlockCache(
    TiledMatrix &matrix,
    size_t block_row_idx,
    size_t block_col_idx)
  : _matrix(matrix) {
    auto& pattern       = matrix.pattern();
    auto  glob_index    = block_row_idx*pattern.blocksize(0)*pattern.extent(1)
                            + block_col_idx*pattern.blocksize(1);
    this->_glob_idx     = glob_index;
    this->_is_local = pattern.is_local(glob_index);
    if (_is_local) {
      _local_ptr = matrix.lbegin() +
                    matrix.pattern().local_index(
                      {block_row_idx*pattern.blocksize(0),
                        block_col_idx*pattern.blocksize(1)}).index;
    }
    this->_size = pattern.blocksize(0) * pattern.blocksize(1);
#ifdef DEBUG
    std::cout << "BlockCache: " << block_row_idx << "x" << block_col_idx << " "
              << " (size=" << _size << ", is_local=" << _is_local
              << ", glob_idx=" << glob_index
              << std::endl;
#endif
  }

  ~LocalBlockCache() {
    if (!this->_is_local && _local_ptr != nullptr) {
#ifdef DEBUG
      std::cout << "Freeing pointer " << _local_ptr << std::endl;
#endif
      free(_local_ptr);
      _local_ptr = nullptr;
    }
  }

  value_t *lbegin() {
    if (!this->_is_local) {
      fetch_data();
    }
    return this->_local_ptr;
  }

  value_t *lend() {
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

  void store() {
    if (!_is_local && this->_local_ptr != NULL) {
      auto begin = this->_matrix.begin() + _glob_idx;
      dash::dart_storage<value_t> ds(_size);
      dart_put_blocking(begin.dart_gptr(), _local_ptr, ds.nelem, ds.dtype);
    }
  }

private:

  void fetch_data() {
    if (this->_local_ptr == nullptr) {
      this->_local_ptr = static_cast<value_t*>(malloc(this->_size * sizeof(value_t)));
      auto begin = this->_matrix.begin() + _glob_idx;
#ifdef DEBUG
      std::cout << "Fetching " << _size << " elements into "
                << _local_ptr << std::endl;
#endif
      dash::dart_storage<value_t> ds(_size);
      dart_get_blocking(
        this->_local_ptr, begin.dart_gptr(), ds.nelem, ds.dtype);
      //dash::copy(begin, end, this->_local_ptr);
    }
  }

private:
  TiledMatrix& _matrix;
  value_t* _local_ptr = nullptr;
  size_t   _size = 0;
  size_t   _glob_idx;
  bool     _is_local  = true;
};


void
compute_single(TiledMatrix& matrix, size_t block_size){

  const size_t num_blocks = matrix.pattern().blockspec().extent(0);
  /**
   * Algorithm taken from
   * https://pm.bsc.es/ompss-docs/examples/01-examples/cholesky
   */

  if (dash::myid() == 0){
    // iterate over column of blocks
    for (size_t k = 0; k < num_blocks; ++k) {

      if (dash::myid() == 0)
        std::cout << "Processing column " << k << std::endl;

      LocalBlockCache block_k(matrix, k, k);

      /**
      * Compute cholesky factorization of block on diagonal
      */
      potrf(block_k.lbegin(), block_size, block_size);
      block_k.store();

      /**
      * Solve the triangular equation system in the block
      */
      for (int i = k+1; i < num_blocks; ++i) {
        LocalBlockCache block_b(matrix, k, i);
        trsm(block_k.lbegin(),
             block_b.lbegin(), block_size, block_size);
        block_b.store();
      }

      // walk to the right
      for (size_t i = k+1; i < num_blocks; ++i) {
        // run down to the diagonal
        LocalBlockCache block_a(matrix, k, i);
        for (size_t j = k+1; j < i; ++j) {
          LocalBlockCache block_c(matrix, j, i);
          LocalBlockCache block_b(matrix, k, j);
          // A[k,i] = A[k,i] - A[k,j] * (A[j,i])^t
          gemm(block_a.lbegin(),
               block_b.lbegin(),
               block_c.lbegin(), block_size, block_size);
          block_c.store();
        }
      }

      // update diagonal blocks
      for (size_t i = k+1; i < num_blocks; ++i) {
        LocalBlockCache block_i(matrix, i, i);
        LocalBlockCache block_a(matrix, k, i);
        // A[j,j] = A[j,j] - A[j,i] * (A[j,i])^t
        syrk(block_a.lbegin(),
             block_i.lbegin(), block_size, block_size);
        block_i.store();
      }
    }
  }

  dash::barrier();
}


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

    LocalBlockCache block_k(matrix, k, k);

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
      LocalBlockCache block_b(matrix, k, i);
      if (block_b.is_local()) {
        trsm(block_k.lbegin(),
             block_b.lbegin(), block_size, block_size);
      }
    }

    dash::barrier();

    // walk to the right
    for (size_t i = k+1; i < num_blocks; ++i) {
      // run down to the diagonal
      LocalBlockCache block_a(matrix, k, i);
      for (size_t j = k+1; j < i; ++j) {
        LocalBlockCache block_c(matrix, j, i);
        if (block_c.is_local()) {
          LocalBlockCache block_b(matrix, k, j);
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
      LocalBlockCache block_i(matrix, i, i);
      if (block_i.is_local()) {
        LocalBlockCache block_a(matrix, k, i);
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

  std::cout << "Matrix: lbegin=" << matrix.lbegin() << ", lend=" << matrix.lend() << std::endl;

  matrix.barrier();
  auto& pattern = matrix.pattern();


  if (dash::myid() == 0) {
    std::cout << "block sizes: " << pattern.blocksize(0) << "x" << pattern.blocksize(1) << std::endl;
    std::cout << "num blocks: "  << pattern.blockspec().extent(0) << "x" << pattern.blockspec().extent(1) << std::endl;
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

#if defined(CHECK_RESULT) && N <= 20
  TiledMatrix matrix_single(N, N, dash::TILE(block_size), dash::TILE(block_size));
  // copy the matrix before compute
  std::copy(matrix.lbegin(), matrix.lend(), matrix_single.lbegin());
  dash::barrier();
  if (dash::myid() == 0)
    print_matrix(matrix_single);
  // compute the correct answer on one unit
  compute_single(matrix_single, block_size);
  if (dash::myid() == 0) {
    std::cout << "########## Expected Result ###############\n";
    print_matrix(matrix_single);
    std::cout << "##########################################\n";
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
void print_matrix(TiledMatrix &matrix)
{
  if (matrix.extent(0) > 100) return;

  std::cout << std::setw(6);

  if (dash::myid() == 0) {
    for (size_t i = 0; i < matrix.extent(0); ++i) {
      for (size_t j = 0; j < matrix.extent(1); ++j) {
        std::cout << std::setw(5) << std::setprecision(3)
                  << static_cast<value_t>(matrix(i, j)) << " ";
      }
      std::cout << std::endl;
    }
  }
}


static
void print_matrix(LocalBlockCache &block, size_t nx, size_t ny)
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


