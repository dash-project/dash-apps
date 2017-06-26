#ifndef DOMAIN_H_INCLUDED
#define DOMAIN_H_INCLUDED

#include <vector>
#include <iostream>
#include <algorithm>

#include <libdash.h>

#include "util.h"
#include "comm.h"
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
  
  NodeMatrixT<Real_t> m_x, m_y, m_z;    // position
  NodeMatrixT<Real_t> m_xd, m_yd, m_zd; // velocity
  NodeMatrixT<Real_t> m_fx, m_fy, m_fz; // force

  NodeMatrixT<Real_t> m_nodalMass;  // mass

  ElemMatrixT<Real_t> m_delv_xi;    // velocity gradient
  ElemMatrixT<Real_t> m_delv_eta;
  ElemMatrixT<Real_t> m_delv_zeta;

  Comm m_comm;
      
 public:
  Domain(std::array<int, 2> nElem,
	 std::array<int, 2> nProcs);

  Comm& comm() { return m_comm; }

  Index_t numProcs() const { return m_comm.numProcs(); } 
  Index_t myRank()   const { return m_comm.myRank(); }

  // number of processes in x and y direction
  Index_t px() const { return m_comm.px(); }
  Index_t py() const { return m_comm.py(); }
  
  Real_t& x(Index_t idx)  { return m_x.lbegin()[idx]; }
  Real_t& y(Index_t idx)  { return m_y.lbegin()[idx]; }
  Real_t& z(Index_t idx)  { return m_z.lbegin()[idx]; }
  Real_t& xd(Index_t idx) { return m_xd.lbegin()[idx]; }
  Real_t& yd(Index_t idx) { return m_yd.lbegin()[idx]; }
  Real_t& zd(Index_t idx) { return m_zd.lbegin()[idx]; }
  Real_t& fx(Index_t idx) { return m_fx.lbegin()[idx]; }
  Real_t& fy(Index_t idx) { return m_fy.lbegin()[idx]; }
  Real_t& fz(Index_t idx) { return m_fz.lbegin()[idx]; }

  Real_t& nodalMass(Index_t idx) { return m_nodalMass.lbegin()[idx]; }
  
  Real_t& delv_xi(Index_t idx)   { return m_delv_xi.lbegin()[idx]; }
  Real_t& delv_eta(Index_t idx)  { return m_delv_eta.lbegin()[idx]; }
  Real_t& delv_zeta(Index_t idx) { return m_delv_zeta.lbegin()[idx]; }

  std::array<Index_t, 2> nNode() const { return m_nNode; }
  std::array<Index_t, 2> nElem() const { return m_nElem; }
  
  Index_t nNode(int dim) const { return m_nNode[dim]; }
  Index_t nElem(int dim) const { return m_nElem[dim]; }

  Index_t maxEdgeSize() const { return std::max(m_nElem[0], m_nElem[1]); }

  // perform data exchange 
  void Exchange(std::array<int,2> dim,
		std::vector<Domain_member>& fields,
		Comm::Scope scope,
		Comm::Direction dir,
		Comm::Action action);

  void PrintNodalMass(int col, int row);
  void PrintForce(int col, int row);
  void PrintPosVel(int col, int row);
  void PrintMonoQ(int col, int row);
};




#endif /* DOMAIN_H_INCLUDED */
