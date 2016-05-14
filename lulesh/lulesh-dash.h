#ifndef LULESH_DASH_H_INCLUDED
#define LULESH_DASH_H_INCLUDED

#include <libdash.h>
#include "lulesh.h"
#include "lulesh-util.h"

/*
 * DASH version of LULESH's 'Domain' data structure.
 *
 * In the MPI version this holds local data in std::vector containers,
 * here we use DASH Matrix. In the MPI version the domain is purely
 * local. Using DASH we have a global view of data from all processes,
 * but of course can also emulate the local view using the .local
 * accessor object.
 */
class Domain
{
private:
  using ElemPatternT = dash::Pattern<3>;
  using NodePatternT = dash::Pattern<3>;

  template<typename ElemType>
  using NodeMatrixT = dash::Matrix<ElemType, 3,
				   NodePatternT::index_type,
				   NodePatternT>;

  template<typename ElemType>
  using ElemMatrixT = dash::Matrix<ElemType, 3,
				   ElemPatternT::index_type,
				   ElemPatternT>;
  // number of local elements and nodes in each dimension
  std::array<Index_t, 3> m_nElem;
  std::array<Index_t, 3> m_nNode;

  // the original code uses the terms 'col', 'row', and 'plane' to
  // refer to the 3D position in the process grid. To get the same
  // ordering, the following equvalencies can be used:
  //
  // m_ts.x() <=> plane
  // m_ts.y() <=> row
  // m_ts.z() <=> col
  dash::TeamSpec<3> m_ts;

  // pattern for element-centered data
  ElemPatternT m_ElemPat;

  // pattern for node-centered data
  NodePatternT m_NodePat;

  NodeMatrixT<Real_t> m_nodalMass;

  NodeMatrixT<Real_t> m_x;  // coordinates
  NodeMatrixT<Real_t> m_y;
  NodeMatrixT<Real_t> m_z;

  ElemMatrixT<Real_t> m_volo ;      // reference volume
  ElemMatrixT<Real_t> m_elemMass;   // mass

  using NodeVecT = std::array<Index_t, 8>;
  // XXX remove me eventually...
  ElemMatrixT<NodeVecT>  m_nodelist ;  // elemToNode connectivity

public:
  Domain(const CmdLineOpts& opts);
  ~Domain();

  void buildMesh();

  Real_t& x(Index_t idx)    { return m_x.lbegin()[idx]; }
  Real_t& y(Index_t idx)    { return m_y.lbegin()[idx]; }
  Real_t& z(Index_t idx)    { return m_z.lbegin()[idx]; }

  // XXX
  Real_t& xd(Index_t idx)   { return m_x.lbegin()[idx]; }
  Real_t& yd(Index_t idx)   { return m_y.lbegin()[idx]; }
  Real_t& zd(Index_t idx)   { return m_z.lbegin()[idx]; }

  // XXX
  Real_t& delv_xi(Index_t idx)   { return m_x.lbegin()[idx]; }
  Real_t& delv_eta(Index_t idx)  { return m_y.lbegin()[idx]; }
  Real_t& delv_zeta(Index_t idx) { return m_z.lbegin()[idx]; }

  Real_t& volo(Index_t idx)         { return      m_volo.lbegin()[idx]; }
  Real_t& elemMass(Index_t idx)     { return  m_elemMass.lbegin()[idx]; }
  Real_t& nodalMass(Index_t idx)    { return m_nodalMass.lbegin()[idx]; }

  Index_t* nodelist(Index_t idx)    { return &m_nodelist.lbegin()[idx][0]; }

  Index_t numElem() const           { return m_nElem[0]*m_nElem[1]*m_nElem[2]; }
  Index_t numElem(size_t dim) const { return m_nElem[dim]; }

  Index_t numNode() const           { return m_nNode[0]*m_nNode[1]*m_nNode[2]; }
  Index_t numNode(size_t dim) const { return m_nNode[dim]; }

  // for compatibility with the old code
  Index_t sizeX() const          { return m_nElem[0]; }
  Index_t sizeY() const          { return m_nElem[1]; }
  Index_t sizeZ() const          { return m_nElem[2]; }

  Index_t colLoc() const         { return m_ts.z(dash::myid()); }
  Index_t rowLoc() const         { return m_ts.y(dash::myid()); }
  Index_t planeLoc() const       { return m_ts.x(dash::myid()); }
  Index_t tp(size_t dim) const   { return m_ts.extent(dim);     }
  Index_t tp() const             { return m_ts.extent(0);       }

  Index_t numRanks() const       { return dash::size(); }
  Index_t maxEdgeSize() const    { return 1+std::max({sizeX(), sizeY(), sizeZ()}); }
  Index_t maxPlaneSize() const   { return maxEdgeSize()*maxEdgeSize(); }

  void print_config(std::ostream& os);
};


Real_t computeChksum(Real_t* ptr, size_t nval);


#endif /* LULESH_DASH_H_INCLUDED */

