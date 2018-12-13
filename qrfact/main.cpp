#include <libdash.h>

#include <iostream>
#include <stdlib.h>
#include <limits>
#include <algorithm>
#include <math.h>
#include <dash/util/Timer.h>

#include <lapacke.h>

#include "common.h"
#include "MatrixBlock.h"

#if defined (DASH_TASKS_COPYIN)
#include "QrFactTasksCopyin.h"
#define USE_TASKS 1
#endif

using value_t = double;
//using PatternT = typename dash::ShiftTilePattern<2>;
using PatternT = typename dash::TilePattern<2>;
//using PatternT = typename dash::TilePattern<2>;
using TiledMatrix = dash::Matrix<value_t, 2, dash::default_index_t, PatternT>;
using Block = MatrixBlock<TiledMatrix>;

//#define DEBUG
//#define CHECK_RESULT
#define FAST_INIT

static
void fill_random(TiledMatrix &matrix);
static
void print_matrix(Block &block, size_t nx, size_t ny);

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

  dash::util::BenchmarkParams bench_params("Cholesky");
  bench_params.print_header();
  bench_params.print_pinning();


#ifdef USE_EXTRAE
  unsigned nvalues = 6;
  Extrae_define_event_type(&et, "Operations", &nvalues, ev, extrae_names);
#endif

  //const size_t  N = 15 * dash::size();

  if (argc < 3) {
    if (dash::myid() == 0) {
      std::cout << argv[0] << " <num_elems> <block_size>\n"
                << "  num_elems:  number of elements in each direction\n"
                << "  block_size: size of block in each direction" << std::endl;
    }
    dash::finalize();
    exit(1);
  }

  size_t num_elem   = atoll(argv[1]); // number of elements in each dimension
  size_t block_size = atoll(argv[2]); // block-size in each dimension


  Timer::Calibrate();

  if (dash::myid() == 0) {
    std::cout << "Allocating matrix: ";
  }

  dash::TeamSpec<2> ts(dash::size(), 1);
  ts.balance_extents();
  TiledMatrix A(num_elem, num_elem,
                     dash::TILE(block_size), dash::TILE(block_size), ts);

  TiledMatrix T(num_elem, num_elem,
                     dash::TILE(block_size), dash::TILE(block_size), ts);

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
  fill_random(A);
#endif
  dash::fill(T.begin(), T.end(), 0.0);

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
              << pattern.blocksize(0) << "x" << pattern.blocksize(1) << std::endl;
    std::cout << "num blocks: "
              << pattern.blockspec().extent(0) << "x"
              << pattern.blockspec().extent(1) << std::endl;
  }

#if defined(CHECK_RESULT)
  decltype(A) A_single(num_elem, num_elem,
                            dash::TILE(block_size), dash::TILE(block_size));
  std::copy(A.lbegin(), A.lend(), A_single.lbegin());
  decltype(T) T_single(num_elem, num_elem,
                            dash::TILE(block_size), dash::TILE(block_size));
  value_t *A_local_copy = new value_t[num_elem*num_elem];
  value_t *A0 = new value_t[num_elem*num_elem];
  // copy the matrix before compute

  dash::barrier();
  if (dash::myid() == 0) {
    std::cout << "########## Initial Matrix ################\n";
    print_matrix(A);
    std::cout << "##########################################\n";
    std::cout << "############ Local Copy ##################\n";
    for (int i = 0; i < num_elem; ++i) {
      for (int j = 0; j < num_elem; ++j) {
        A_local_copy[i*num_elem + j] = A0[i*num_elem + j] = A_single(i, j);
      }
    }
    print_matrix(A_local_copy, num_elem);

    std::cout << "##########################################\n";
    int num_elem_i = num_elem;
    int lwork = num_elem*num_elem;
    int info;
    value_t *tau  = new value_t[lwork];
    //LAPACK_dgeqrf(&num_elem_i, &num_elem_i, A_local_copy, &num_elem_i, tau, work, &lwork, &info);
    int ret = LAPACKE_dgeqrt(LAPACK_ROW_MAJOR, num_elem_i, num_elem_i, block_size, A_local_copy, num_elem, tau, num_elem);
    delete[] tau;
    std::cout << "LAPACKE_dgeqrt returned " << ret << std::endl;
    for (int i = 0; i < num_elem; ++i) {
      for (int j = 0; j < num_elem; ++j) {
        A_single(i, j) = A_local_copy[i*num_elem + j];
      }
    }
  }
  dash::barrier();
#if 0
  matrix.barrier();
  if (dash::myid() == 0 && num_elem <=20) {
    if (num_elem <=20) {
      print_matrix(A_single);
    }
    std::cout << "Computing verification matrix" << std::endl;
  }
  // compute the correct answer on one unit
  compute_single(A_single, T_single, block_size);
#endif
  if (dash::myid() == 0) {
    std::cout << "Done computing verification matrix" << std::endl;
    if (num_elem <=20) {
      std::cout << "########## Expected Result ###############\n";
      print_matrix(A_single);
      std::cout << "##########################################\n";
    }
  }
#endif

  Timer t;
  compute(A, T, block_size);
  //compute_single(A, T, block_size);
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
  verify_matrix(A, A_single);
  value_t *T_local_copy = new value_t[num_elem*num_elem];
  value_t *B  = new value_t[num_elem*num_elem];
  value_t *X  = new value_t[num_elem*num_elem];
  value_t *Id = new value_t[num_elem*num_elem];
  value_t *work = new value_t[num_elem];
  double eps = LAPACKE_dlamch_work('e');
  for (int i  = 0; i < num_elem; ++i) {
    for (int j = 0; j < num_elem; ++j) {
      A_local_copy[i*num_elem + j] = A(i, j);
      T_local_copy[i*num_elem + j] = T(i, j);
    }
  }
  int k = num_elem / block_size;
  value_t *Q = new value_t[num_elem*num_elem];
  CORE_dlacpy(PlasmaUpper,      num_elem, num_elem, A_local_copy, num_elem, Q, num_elem);
  LAPACKE_dorgqr(LAPACK_ROW_MAJOR, num_elem, num_elem, k, Q, num_elem, T_local_copy);
  CORE_dplrnt(num_elem, num_elem, X, num_elem, num_elem, 0, 0, 2354);
  CORE_dlacpy(PlasmaUpperLower, num_elem, num_elem, X, num_elem, B, num_elem);
  // check diagonality
  CORE_dlaset(PlasmaUpperLower, num_elem, num_elem, 0., 1., Id, num_elem);
  CORE_dsyrk(PlasmaUpper, PlasmaTrans, num_elem, num_elem, 1.0, Q, num_elem, -1.0, Id, num_elem);
  value_t norm;
  CORE_dlansy(PlasmaInfNorm, PlasmaUpper, num_elem, Id, num_elem, work, &norm);
  value_t result = norm / (num_elem / eps);
  if (isnan(result) || isinf(result) || result > 60.0) {
    std::cout << "Norm looks suspicious: " << result << std::endl;
  } else {
    std::cout << "Diagonality is correct!" << std::endl;
  }

  // check factorization
  value_t *Residual = new value_t[num_elem*num_elem];
  value_t *R        = new value_t[num_elem*num_elem];
  CORE_dlacpy(PlasmaUpperLower, num_elem, num_elem, A0, num_elem, Residual, num_elem);
  CORE_dlaset(PlasmaUpperLower, num_elem, num_elem, 0., 0., R, num_elem);
  CORE_dlacpy(PlasmaUpper,      num_elem, num_elem, A_local_copy, num_elem, R, num_elem);
  /* Compute Residual = A0 - Q*R */
  CORE_dgemm(PlasmaNoTrans, PlasmaNoTrans, num_elem, num_elem, num_elem, -1.0, Q, num_elem, R, num_elem, 1.0, Residual, num_elem);
  value_t Anorm, Rnorm;
  CORE_dlange(PlasmaInfNorm, num_elem, num_elem, Residual, num_elem, work, &Rnorm);
  CORE_dlange(PlasmaInfNorm, num_elem, num_elem, A0, num_elem, work, &Anorm);
  result = Rnorm / ( Anorm * num_elem * eps);

  if (isnan(result) || isinf(result) || result > 60.0) {
    std::cout << "Norm looks suspicious: " << result << std::endl;
  } else {
    std::cout << "Factorization is correct!" << std::endl;
  }

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

