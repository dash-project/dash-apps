#include <iostream>
#include <libdash.h>
#include "dash/Init.h"
#include "dash/halo/Halo.h"
#include "dash/halo/Types.h"

#include "minimon.h"

using namespace std;

constexpr dash::dim_t DIMENSION = 2;
constexpr size_t INIT_WIDTH = 5;

using PatternT    = dash::Pattern<DIMENSION>;
using MatrixT     = dash::Matrix<double, DIMENSION, typename PatternT::index_type, PatternT>;
using StencilT     = dash::halo::StencilPoint<DIMENSION>;
using StencilSpecT = dash::halo::StencilSpec<StencilT,4>;
using GlobBoundSpecT     = dash::halo::GlobalBoundarySpec<DIMENSION>;
using HaloMatrixWrapperT = dash::halo::HaloMatrixWrapper<MatrixT>;

using array_t      = dash::Array<double>;

MiniMon minimon{};

constexpr double dx = 1.0;
constexpr double dy = 1.0;
constexpr double dt = 0.05;
constexpr double k = 1.0;

double calcEnergy(const MatrixT& m, array_t &a) {
  *a.lbegin() = std::accumulate( m.lbegin(), m.lend(), 0.0);
  a.barrier();
  double energy = 0.0;
  if(dash::myid() == 0)
    energy = std::accumulate(a.begin(), a.end(), 0.0);
  a.barrier();

  return energy;
}

int main(int argc, char *argv[])
{
  minimon.enter();

  if (argc < 3) {
    cerr << "Not enough arguments ./<prog> matrix_ext iterations" << endl;
    return 1;
  }

  auto matrix_ext = std::atoi(argv[1]);
  auto iterations = std::atoi(argv[2]);

  minimon.enter();
  dash::init(&argc, &argv);
  minimon.leave("initialize");

  auto myid = dash::myid();
  auto ranks = dash::size();

  dash::DistributionSpec<DIMENSION> dist(dash::BLOCKED, dash::BLOCKED);
  dash::TeamSpec<DIMENSION> tspec;
  tspec.balance_extents();
  PatternT pattern(dash::SizeSpec<DIMENSION>(matrix_ext, matrix_ext),dist, tspec, dash::Team::All());
  auto src_matrix = MatrixT(pattern);
  auto dst_matrix = MatrixT(pattern);

  StencilSpecT stencil_spec( StencilT(-1, 0), StencilT(1, 0),
                             StencilT( 0,-1), StencilT(0, 1));

  // Periodic/cyclic global boundary values for both dimensions
  GlobBoundSpecT bound_spec(dash::halo::BoundaryProp::CUSTOM,dash::halo::BoundaryProp::CUSTOM);

  // HaloWrapper for source and destination partitions
  HaloMatrixWrapperT src_halo(src_matrix, bound_spec, stencil_spec);
  HaloMatrixWrapperT dst_halo(dst_matrix, bound_spec, stencil_spec);
  auto src_halo_ptr = &src_halo;
  auto dst_halo_ptr = &dst_halo;

  // Stencil specific operator for both partitions
  auto src_stencil_op = src_halo_ptr->stencil_operator(stencil_spec);
  auto dst_stencil_op = dst_halo_ptr->stencil_operator(stencil_spec);
  auto src_op_ptr = &src_stencil_op;
  auto dst_op_ptr = &dst_stencil_op;

  for(auto i = 0; i < src_matrix.local.extent(0); ++i) {
    for(auto j = 0; j < src_matrix.local.extent(1); ++j) {
      if(i < INIT_WIDTH  && j < INIT_WIDTH) {
        src_matrix.local[i][j] = 1;
        dst_matrix.local[i][j] = 1;
      } else {
        src_matrix.local[i][j] = 0;
        dst_matrix.local[i][j] = 0;
      }
    }
  }

  src_halo.set_custom_halos([](const auto& coords) { return 0;});
  dst_halo.set_custom_halos([](const auto& coords) { return 0;});


  // initial total energy
  auto energy = array_t(ranks);
  double initEnergy = calcEnergy(src_halo_ptr->matrix(), energy);

  src_matrix.barrier();
  minimon.enter();

  for (auto i = 0; i < iterations; ++i) {
    minimon.enter();
    auto dst_matrix_lbegin = dst_halo_ptr->matrix().lbegin();

    minimon.enter();
    src_halo_ptr->update_async(); // start asynchronous Halo data exchange
    minimon.leave("async");

    minimon.enter();
    // Calculation of all inner partition elements
    auto iend = src_op_ptr->inner.end();
    for(auto it = src_op_ptr->inner.begin(); it != iend; ++it) {
      auto center = *it;
      double dtheta = (it.value_at(0) + it.value_at(1) - 2 * center) / (dx * dx) +
                      (it.value_at(2) + it.value_at(3) - 2 * center) / (dy * dy);
      dst_matrix_lbegin[it.lpos()] = center + k * dtheta * dt;
    }
    minimon.leave("inner");

    minimon.enter();
    src_halo_ptr->wait(); // Wait until all Halo data exchanges are finished
    minimon.leave("wait");

    minimon.enter();
    // Calculation of all boundary partition elements
    auto bend = src_op_ptr->boundary.end();
    for(auto it = src_op_ptr->boundary.begin(); it != bend; ++it) {
      auto center = *it;
      double dtheta = (it.value_at(0) + it.value_at(1) - 2 * center) / (dx * dx) +
                      (it.value_at(2) + it.value_at(3) - 2 * center) / (dy * dy);
      dst_matrix_lbegin[it.lpos()] = center + k * dtheta * dt;
    }
    minimon.leave("boundary");

    // swap source and destination partitions and operators
    std::swap(src_halo_ptr, dst_halo_ptr);
    std::swap(src_op_ptr, dst_op_ptr);
    minimon.leave("calc iter");
  }

  dash::barrier();
  minimon.leave("calc total");

  // final total energy
  double endEnergy = calcEnergy(src_halo_ptr->matrix(), energy);

  // Output
  if (myid == 0) {
    cout << fixed;
    cout.precision(5);
    cout << "InitEnergy=" << initEnergy << endl;
    cout << "EndEnergy=" << endEnergy << endl;
    cout << "DiffEnergy=" << endEnergy - initEnergy  << endl;
    cout << "Matrixspec: " << matrix_ext << " x " << matrix_ext << endl;
    cout << "Iterations: " << iterations << endl;
    cout.flush();
  }

  dash::Team::All().barrier();
  minimon.leave("total");

  if(myid == 0)
    std::cout << "# unit_id;function_name;num_calls;avg_runtime;min_runtime;max_runtime"
              << std::endl;

  for(auto i = 0; i < dash::size(); ++i) {
    if(i == myid)
      minimon.print(myid);
    dash::Team::All().barrier();
  }

  dash::finalize();

  return 0;
}
