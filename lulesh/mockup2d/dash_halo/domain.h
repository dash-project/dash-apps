#ifndef DOMAIN_H_INCLUDED
#define DOMAIN_H_INCLUDED

#include <vector>
#include <iostream>
#include <algorithm>

#include <libdash.h>
#include <dash/experimental/HaloMatrixWrapper.h>
#include "util.h"
#include "lulesh-mock.h"


/*

 Mock-up of the communication behavior of LULESH (2D version)
 ============================================================

      /  z, planes (not included in 2D version)
     /
    +-----> x, cols
    |
    |
    v y, rows

  Each process/unit is responsible for a part of the overall
  simulation, called 'domain'. The size of the local domain is given
  by dx (number of columns) and dy (number of rows):

       |---------- dx -----------|

    -  A#########################B <- Y0 (top row)
    |  # --- storage order --->  #
    |  #                         #
  dy|  #                         #
    |  #                         #
    |  #                         #
    -  C#########################D <- Y1 (bottom row)
       ^                         ^
       X0 (left col)             X1 (right col)

       A = X0Y0: top-left corner
       B = X1Y0: top-right corner
       C = X0Y1: bottom-left corner
       D = X1Y1: bottom-right corner

  Note that the storage order in the above sketch is row-by-row but
  that rows correspond to the first (x) dimension, so in fact this
  corresponds to a column-major 2D array/matrix. This was chosen such
  that the transition to 3D is easier and more natural.

  Edges and corners (+planes in 3D) are communicated with neighboring
  processes. To simplify the implemenation compared to the original
  implementation 'tags' as shown above are used to denote what is
  communicated.

  X0, X1, Y0, Y1         denote edges
  X0Y0, X0Y1, X1Y0, X1Y1 denote corners

*/

class Domain
{
 private:
  using ElemPatternT = dash::Pattern<2, dash::COL_MAJOR>;
  using NodePatternT = dash::Pattern<2, dash::COL_MAJOR>;
  using StencilT = dash::experimental::Stencil<2>;
  using StencilSpecT = dash::experimental::StencilSpec<2, 8>;
  using CycleSpecT = dash::experimental::CycleSpec<2>;
  using Cycle     = dash::experimental::Cycle;

  template<typename ElemType>
  using NodeMatrixT = dash::Matrix<ElemType, 2,
				   NodePatternT::index_type,
				   NodePatternT>;

  template<typename ElemType>
  using ElemMatrixT = dash::Matrix<ElemType, 2,
				   ElemPatternT::index_type,
				   ElemPatternT>;

  // number of local elements and nodes in each dim.
  std::array<Index_t, 2> m_nElem;
  std::array<Index_t, 2> m_nNode;

  dash::TeamSpec<2> m_teamspec;

  // pattern for element-centered data
  ElemPatternT m_ElemPat;

  // pattern for node-centered data
  NodePatternT m_NodePat;

  NodeMatrixT<TripleT> m_pos;    // position
  NodeMatrixT<TripleT> m_vel; // velocity
  NodeMatrixT<TripleT> m_force; // force

  NodeMatrixT<Real_t> m_nodalMass;  // mass

  ElemMatrixT<TripleT> m_delv;    // velocity gradient

  StencilSpecT stencil_spec;

  dash::experimental::HaloMatrixWrapper<NodeMatrixT<TripleT>,StencilSpecT> halo_pos;
  dash::experimental::HaloMatrixWrapper<NodeMatrixT<TripleT>,StencilSpecT> halo_vel;
  dash::experimental::HaloMatrixWrapper<NodeMatrixT<TripleT>,StencilSpecT> halo_force;
  dash::experimental::HaloMatrixWrapper<NodeMatrixT<Real_t>,StencilSpecT>  halo_nodalMass;
  dash::experimental::HaloMatrixWrapper<ElemMatrixT<TripleT>,StencilSpecT> halo_delv;
  //Comm m_comm;
  std::array<int, 2> _nProcs;

 public:
  Domain(const std::array<int, 2>& nElem,
	 const std::array<int, 2>& nProcs);

  Index_t numProcs() const { return dash::size(); }
  Index_t myRank()   const { return dash::myid(); }

  // number of processes in x and y direction
  Index_t px() const { return _nProcs[0]; }
  Index_t py() const { return _nProcs[1]; }


   // WORKAROUND because DASH currently does not
   // allow COL_MAJOR teamspec:
   // row/col and px() py() are switched
  Index_t col() const { return myRank()/py(); }
  Index_t row() const { return myRank()%py(); }

  Real_t& x(Index_t idx)  { return m_pos.lbegin()[idx].x; }
  Real_t& y(Index_t idx)  { return m_pos.lbegin()[idx].y; }
  Real_t& z(Index_t idx)  { return m_pos.lbegin()[idx].z; }
  Real_t& xd(Index_t idx) { return m_vel.lbegin()[idx].x; }
  Real_t& yd(Index_t idx) { return m_vel.lbegin()[idx].y; }
  Real_t& zd(Index_t idx) { return m_vel.lbegin()[idx].z; }
  Real_t& fx(Index_t idx) { return m_force.lbegin()[idx].x; }
  Real_t& fy(Index_t idx) { return m_force.lbegin()[idx].y; }
  Real_t& fz(Index_t idx) { return m_force.lbegin()[idx].z; }

  Real_t& nodalMass(Index_t idx) { return m_nodalMass.lbegin()[idx]; }

  Real_t& delv_xi(Index_t idx)   { return m_delv.lbegin()[idx].x; }
  Real_t& delv_eta(Index_t idx)  { return m_delv.lbegin()[idx].y; }
  Real_t& delv_zeta(Index_t idx) { return m_delv.lbegin()[idx].z; }

  std::array<Index_t, 2> nNode() const { return m_nNode; }
  std::array<Index_t, 2> nElem() const { return m_nElem; }

  Index_t nNode(int dim) const { return m_nNode[dim]; }
  Index_t nElem(int dim) const { return m_nElem[dim]; }

  Index_t maxEdgeSize() const { return std::max(m_nElem[0], m_nElem[1]); }

  // perform data exchange
  void ExchangePos();
  void ExchangeVelocity();
  void exchangeForce();
  void exchangeNodalMass();
  void ExchangeGradiant();

  void PrintNodalMass(int col, int row);
  void PrintForce(int col, int row);
  void PrintPosVel(int col, int row);
  void PrintMonoQ(int col, int row);

 private:
  template<typename HaloMatrixWrapperT>
  void exchangeAdd(HaloMatrixWrapperT halo_wrapper);
};




#endif /* DOMAIN_H_INCLUDED */
