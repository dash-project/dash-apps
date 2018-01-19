#include <libdash.h>

#include <iostream>
#include <stdlib.h>
#include <limits>
#include <algorithm>
#include <dash/util/Timer.h>

#include "common.h"
#include "MatrixBlock.h"

#if defined (DASH_TASKS)
#include "CholeskyTasks.h"
#define USE_TASKS 1
#elif defined (DASH_TASKS_PREFETCH)
#include "CholeskyTasksPrefetch.h"
#define USE_TASKS 1
#elif defined (DASH_TASKS_COPYIN)
#include "CholeskyTasksCopyIn.h"
#define USE_TASKS 1
#elif defined (DASH_PREFETCH)
#include "CholeskyPrefetch.h"
#define USE_TASKS 0
#elif defined (DASH_OMP_PREFETCH)
#include "CholeskyOMPPrefetch.h"
#define USE_TASKS 0
#elif defined (MPI_OMP)
#include "CholeskyMPIOMP.h"
#define USE_TASKS 0
#else
#include "Cholesky.h"
#define USE_TASKS 0
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
  TiledMatrix matrix(num_elem, num_elem,
                     dash::TILE(block_size), dash::TILE(block_size), ts);

  /*
  TiledMatrix matrix(num_elem, num_elem,
                     dash::NONE, dash::BLOCKCYCLIC(block_size), dash::TeamSpec<2>(1, dash::size()));
  */
  if (dash::myid() == 0) {
    std::cout << "Done." << std::endl;
  }
  fill_random(matrix);

  matrix.barrier();
  auto& pattern = matrix.pattern();


  if (dash::myid() == 0) {
    std::cout << "Implementation: " << CHOLESKY_IMPL << std::endl;
    std::cout << "Matrix: " << num_elem << "x" << num_elem << std::endl;
    std::cout << "block sizes: "
              << pattern.blocksize(0) << "x" << pattern.blocksize(1) << std::endl;
    std::cout << "num blocks: "
              << pattern.blockspec().extent(0) << "x"
              << pattern.blockspec().extent(1) << std::endl;
  }

#if defined(CHECK_RESULT)
  decltype(matrix) matrix_single(num_elem, num_elem,
                            dash::TILE(block_size), dash::TILE(block_size));
  // copy the matrix before compute
  std::copy(matrix.lbegin(), matrix.lend(), matrix_single.lbegin());
  matrix.barrier();
  if (dash::myid() == 0 && num_elem <=20) {
    if (num_elem <=20) {
      print_matrix(matrix_single);
    }
    std::cout << "Computing verification matrix" << std::endl;
  }
  // compute the correct answer on one unit
  compute_single(matrix_single, block_size);
  if (dash::myid() == 0) {
    std::cout << "Done computing verification matrix" << std::endl;
    if (num_elem <=20) {
      std::cout << "########## Expected Result ###############\n";
      print_matrix(matrix_single);
      std::cout << "##########################################\n";
    }
  }
#endif

  Timer t;
  compute(matrix, block_size);
  dash::barrier();
  auto elapsed = t.Elapsed(); // time in us
  double flops = (1.0 / 3.0) * num_elem * num_elem * num_elem;

  if (dash::myid() == 0) {
    std::cout << "Cholesky factorization of " << num_elem << "x" << num_elem
              << " done after " << elapsed/1E3 << "ms ("
              << flops/elapsed/1E3 << " GF/s)" << std::endl;
  }

  if (num_elem <= 20)
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
#elif defined(FAST_INIT) && USE_TASKS
  using value_t = typename TiledMatrix::value_type;
  dash::tasks::parallel_for(matrix.lbegin(), matrix.lend(),
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
        dlarnv_(&intONE, &ISEED[0], &n, lbegin + num_elem);
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


