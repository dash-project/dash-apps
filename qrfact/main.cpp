#include <libdash.h>

#include <iostream>
#include <stdlib.h>
#include <limits>
#include <algorithm>
#include <math.h>
#include <dash/util/Timer.h>

//#include <lapacke.h>

#include <mkl.h>

#include "common.h"
#include "MatrixBlock.h"

#if defined (DASH_TASKS_COPYIN)
#include "QrFactTasksCopyin.h"
#define USE_TASKS 1
#endif

using value_t = double;
//using PatternT = typename dash::ShiftTilePattern<2>;
using PatternT = typename dash::TilePattern<4>;
//using PatternT = typename dash::TilePattern<2>;
using TiledMatrix = dash::Matrix<value_t, 4, dash::default_index_t, PatternT>;
using Block = MatrixBlock<TiledMatrix>;

//#define DEBUG
//#define CHECK_RESULT
#define FAST_INIT
//#define CHECK_RESULT_GEQRT

static
void fill_random(TiledMatrix &matrix);
static
void print_matrix(Block &block, size_t nx, size_t ny);

template<typename ValueT>
static int check_orthogonality(ValueT *matrix, int num_elem);

template<typename ValueT>
static int check_factorization(ValueT *A0, ValueT *A, ValueT *Q, size_t N);

int main(int argc, char **argv)
{
  dash::init(&argc, &argv);

  if (!dash::is_multithreaded()) {
    if (dash::myid() == 0) {
      std::cerr << "ERROR: Missing support for multi-threading!" << std::endl;
    }
    dash::finalize();
    exit(1);
  }

  dash::util::BenchmarkParams bench_params("QR Factorization");
  bench_params.print_header();
  bench_params.print_pinning();

  if (argc < 3) {
    if (dash::myid() == 0) {
      std::cout << argv[0] << " <num_elems> <block_size> [<sbm = 4> [<sbn=1>]]\n"
                << "  num_elems:  number of elements in each direction\n"
                << "  block_size: size of block in each direction" 
                << "  sbm: super-block size in M"
                << "  sbn: super-block size in N" << std::endl;
    }
    dash::finalize();
    exit(1);
  }

  size_t num_elem   = atoll(argv[1]); // number of elements in each dimension
  size_t block_size = atoll(argv[2]); // block-size in each dimension
  size_t sb_size_m = 4;
  size_t sb_size_n = 1;
  if (argc > 3) {
    sb_size_m = sb_size_n = atoll(argv[3]);
  }
  if (argc > 4) {
    sb_size_n = atoll(argv[4]);
  }


  Timer::Calibrate();

  if (dash::myid() == 0) {
    std::cout << "Allocating matrix: \n";
  }

  dash::TeamSpec<2> ts2d(dash::size(), 1);
  ts2d.balance_extents();
  dash::TeamSpec<4> ts(ts2d.extent(0), ts2d.extent(1), 1, 1);
  size_t num_blocks = num_elem / block_size;
  std::cout << "num_blocks: " << num_blocks << ", SB " << sb_size_m << "x" << sb_size_n << ", block_size " << block_size << std::endl;
  TiledMatrix A(num_blocks, num_blocks,         // super-blocks
                dash::TILE(sb_size_m), dash::TILE(sb_size_n),
                block_size, block_size,
                dash::TILE(block_size), dash::TILE(block_size), ts);

  TiledMatrix T(num_blocks, num_blocks,         // super-blocks
                dash::TILE(sb_size_m), dash::TILE(sb_size_n),
                block_size, block_size,
                dash::TILE(block_size), dash::TILE(block_size), ts);

  std::cout << dash::myid() << ": extents " << A.local_size() << std::endl;

#if defined(CHECK_RESULT)
#if 0
  for (auto iter = A.lbegin(); iter != A.lend(); ++iter) {
    *iter = dash::myid() + c++ / 1E4;
  }
#else
  int block = 0;
  for (size_t i = 0; i < num_blocks; ++i) {
    for (size_t j = 0; j < num_blocks; ++j) {
      if (A(i, j, 0, 0).is_local()) {
        int elem = 0;
        for (size_t k = 0; k < block_size; ++k) {
          for (size_t l = 0; l < block_size; ++l) {
            A(i, j, k, l) = dash::myid() + block / 1E2 + elem++ / 1E4;
          }
        }
        ++block;
      }
    }
  }
  dash::barrier();
  print_matrix(A);
  dash::barrier();
#endif

  value_t *A0 = new value_t[num_elem*num_elem];
  {
  size_t c = 0;
  for (size_t i = 0; i < num_blocks; ++i) {
    for (size_t k = 0; k < block_size; ++k) {
      for (size_t j = 0; j < num_blocks; ++j) {
        for (size_t l = 0; l < block_size; ++l) {
          A0[c++] = A(i, j, k, l);
        }
      }
    }
  }
  }

  print_matrix(A0, num_elem);
  dash::barrier();
#endif // defined(CHECK_RESULT)



#ifdef CHECK_RESULT_GEQRT

  if (dash::myid() == 0) {
    value_t *A_local = new value_t[num_elem*num_elem];
    value_t *T_local = new value_t[num_elem*num_elem];

    value_t *work    = new double[block_size*block_size];
    value_t *tau     = new double[block_size];

    for (size_t i = 0; i < num_blocks; ++i) {
      for (size_t j = 0; j < num_blocks; ++j) {
        for (size_t k = 0; k < block_size; ++k) {
          for (size_t l = 0; l < block_size; ++l) {
            A_local[i*num_blocks*block_size*block_size + j*block_size*block_size + k*block_size + l] = A(i, j, k, l);
          }
        }
      }
    }
    std::cout << "LOCAL VERIFICATION MATRIX A_LOCAL" << std::endl;
    geqrt(A_local, T_local, tau, work, block_size);
    std::cout << "LOCAL VERIFICATION MATRIX A_LOCAL AFTER GEQRT" << std::endl;
    print_matrix(A_local, num_elem);
    delete[] A_local;
    delete[] T_local;
    delete[] work;
    delete[] tau;
  }


#endif // CHECK_RESULT_GEQRT  

#if 0
  dash::barrier();
  print_matrix(A);
  std::cout << std::endl;
  for (size_t i = 0; i < num_blocks; ++i) {
    for (size_t j = 0; j < num_blocks; ++j) {
      std::cout << "Block " << i << ", " << j << ": " << std::endl;
      if (A(i, j, 0, 0).is_local()) {
        print_block(A(i, j, 0, 0).local_ptr(), block_size);
      }
    }
  }
  dash::barrier();
#endif

  /*
  TiledMatrix matrix(num_elem, num_elem,
                     dash::NONE, dash::BLOCKCYCLIC(block_size), dash::TeamSpec<2>(1, dash::size()));
  */
  if (dash::myid() == 0) {
    std::cout << "Done." << std::endl;
  }

//  if (dash::size() == 1)
#ifdef CHECK_RESULT
  {
    int num_blocks = num_elem / block_size;
    std::cout << "Generating matrix using DPLRNT" << std::endl;
    for (int i = 0; i < num_blocks; ++i) {
      for (int j = 0; j < num_blocks; ++j) {
        Block b(A, i, j);
        if (b.is_local())
          CORE_dplrnt(block_size, block_size, b.lbegin(), block_size, num_elem, block_size*i, block_size*j, 3872);
      }
    }
  }
#else
  //fill_random(A);
  dash::fill(A.begin(), A.end(), 0.0);
#endif
  //dash::fill(T.begin(), T.end(), 0.0);

#if 0
    using LocalBlockCache = MatrixBlock<TiledMatrix>;

    for (int i = 0; i < num_elem/block_size; ++i) {
      for (int j = 0; j < num_elem/block_size; ++j) {
        LocalBlockCache block(matrix, i, j);
        if (block.is_local()) {
          std::fill(block.lbegin(), block.lend(), (double)(dash::myid()*10 + i + ((double)j)/10));
        }
      }
    }
#endif

  A.barrier();
  auto& pattern = A.pattern();

  if (dash::myid() == 0) {
    std::cout << "Implementation: " << QRFACT_IMPL << std::endl;
    std::cout << "Matrix: " << num_elem << "x" << num_elem << std::endl;
    std::cout << "block sizes: "
              << block_size << "x" << block_size << std::endl;
    std::cout << "num blocks: "
              << pattern.extent(0) << "x"
              << pattern.extent(1) << std::endl;
  }

  Timer t;
  compute(A, T, block_size);
  dash::barrier();
  auto elapsed = t.Elapsed(); // time in us
  double flops = (4.0 / 3.0) * num_elem * num_elem * num_elem;


  if (dash::myid() == 0) {
    std::cout << "QR factorization of " << num_elem << "x" << num_elem
              << " done after " << elapsed/1E3 << "ms ("
              << flops/elapsed/1E3 << " GF/s)" << std::endl;
  }

  if (num_elem <= 20)
    print_matrix(A);

#if defined(CHECK_RESULT)
  value_t *A_local_copy = new value_t[num_elem*num_elem];
  value_t *T_local_copy = new value_t[num_elem*num_elem];
  value_t *B            = new value_t[num_elem*num_elem];
  value_t *X            = new value_t[num_elem*num_elem];
  double eps = LAPACKE_dlamch_work('e');
  size_t c = 0;
  for (size_t i = 0; i < num_blocks; ++i) {
    for (size_t k = 0; k < block_size; ++k) {
      for (size_t j = 0; j < num_blocks; ++j) {
        for (size_t l = 0; l < block_size; ++l) {
          A_local_copy[c] =  A(i, j, k, l);
          T_local_copy[c] =  T(i, j, k, l);
          c++;
        }
      }
    }
  }
  //geqrt(A_local, T_local, tau, work, block_size);
  int k = num_elem / block_size;
  size_t N = num_elem;
  value_t *Q = new value_t[N*N];
  CORE_dlacpy(PlasmaUpperLower,    N, N, A_local_copy, N, Q, N);
  LAPACKE_dorgqr(LAPACK_ROW_MAJOR, N, N, N, Q, N, T_local_copy);
  CORE_dplrnt(N, N, X, N, N, 0, 0, 2354);
  CORE_dlacpy(PlasmaUpperLower, N, N, X, N, B, N);
  // check diagonality
  check_orthogonality(Q, N);
  // check factorization
  check_factorization(A0, A_local_copy, Q, N);

  // clean up
  delete[] Q;
#endif

  dash::finalize();

  return 0;
}

static
void fill_random(TiledMatrix &matrix)
{
#if defined(DEBUG) //|| defined(CHECK_RESULT)
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
#elif defined(FAST_INIT) && USE_TASKS
  using value_t = typename TiledMatrix::value_type;
  dash::tasks::taskloop(matrix.lbegin(), matrix.lend(),
    [](value_t* first, value_t* last){
      int ISEED[4] = {0,0,0,1};
      int intONE=1;
      size_t num_elem_total = last - first;
      size_t num_elem = 0;
      while (num_elem < num_elem_total) {
        size_t to_copy = num_elem_total - num_elem;
        int n = std::min(to_copy,
                         static_cast<size_t>(std::numeric_limits<int>::max()));
        dlarnv_(&intONE, &ISEED[0], &n, first + num_elem);
        num_elem += n;
      }
    });
  dash::tasks::complete();
#elif defined(FAST_INIT) && defined(_OPENMP)
  constexpr int rand_max = 100;
#pragma omp parallel
#pragma omp single
{
  size_t pos = 0;
  size_t chunk_size = matrix.local_size() / omp_get_num_threads();
  size_t num_chunks = matrix.local_size() / chunk_size;
  for (pos = 0; pos < matrix.local_size(); pos += chunk_size) {
#pragma omp task
{
      int ISEED[4] = {0,0,0,1};
      int intONE=1;
      auto lbegin = matrix.lbegin() + pos;
      size_t num_elem_total = std::min(chunk_size, matrix.local_size() - pos);
      size_t num_elem = 0;
      while (num_elem < num_elem_total) {
        size_t to_copy = num_elem_total - num_elem;
        int n = std::min(to_copy,
                         static_cast<size_t>(std::numeric_limits<int>::max()));
        LAPACK_dlarnv(intONE, ISEED, n, lbegin + num_elem);
        num_elem += n;
      }
}
  }
}
#else
  constexpr int rand_max = 100;
  for (auto it = matrix.lbegin(); it != matrix.lend(); ++it) {
    *it = (rand())%(rand_max);
  }
#endif
  // add to diagonal
  size_t ext = matrix.extent(0);
  for (size_t i = 0; i < ext; ++i) {
    if (matrix(i, i).is_local()) matrix(i, i) = ext;
  }
}

template<typename ValueT>
static int check_orthogonality(ValueT *A, int N)
{
  size_t NN = N*N;
  double eps = LAPACKE_dlamch_work('e');

  // create identity matrix
  ValueT *Id = new ValueT[NN];
  std::fill(Id, Id+NN, 0.0);
  for (int i = 0; i < NN; ++i) {
    for (int j = 0; j < NN; ++j) {
      if (i == j) Id[i*N+j] = 1.0;
    }
  }

  // compute I - Q'Q
  cblas_dsyrk(CblasRowMajor, CblasUpper, CblasTrans, N, N, 1.0, A, N, -1.0, Id, N);

  ValueT *work = new ValueT[N];
  // compute infinity norm
  ValueT normQ;
  CORE_dlansy('i', 'u', N, Id, N, work, &normQ);
  delete[] work;
  delete[] Id;

  ValueT result = normQ / (N * eps);

  if ( isnan(result) || isinf(result) || (result > 60.0) ) {
     std::cout << "-- Orthogonality is suspicious !" << std::endl;
  }
  else {
    std::cout << "-- Orthogonality is CORRECT !" << std::endl;
  }

  return 0;
}

template<typename ValueT>
static int check_factorization(ValueT *A0, ValueT *A, ValueT *Q, size_t N)
{
  size_t NN = N*N;
  ValueT eps = LAPACKE_dlamch_work('e');
  ValueT *R        = new ValueT[NN];
  ValueT *Residual = new ValueT[NN];

  // copy original A into residual
  CORE_dlacpy(PlasmaUpperLower, N, N, A0, N, Residual, N);

  // extract R
  std::fill(R, R+NN, 0.0);
  CORE_dlacpy(PlasmaUpper, N, N, A, N, R, N);
  
  // perform Residual = Aorig - Q*R
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, -1.0, Q, N, R, N, 1.0, Residual, N);

  // compute Residual norm
  ValueT normR;
  ValueT *work = new ValueT[N];
  CORE_dlange(PlasmaInfNorm, N, N, Residual, N, work, &normR);

  // compute Aorig norm
  ValueT normA;
  CORE_dlange(PlasmaInfNorm, N, N, A0, N, work, &normA);

  ValueT result = normR / ( normA * N * eps);

  if ( isnan(result) || isinf(result) || (result > 60.0) ) {
     std::cout << "-- Factorization is suspicious!" << std::endl;
  } else { 
     std::cout << "-- Factorization is CORRECT !" << std::endl;
  }

  delete[] work;
  delete[] R;
  delete[] Residual;
  return 0;
}

