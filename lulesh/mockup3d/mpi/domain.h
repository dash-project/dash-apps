#ifndef DOMAIN_H_INCLUDED
#define DOMAIN_H_INCLUDED

#include <vector>
#include <algorithm>

#include "comm.h"
#include "lulesh-mock.h"


/*

 Mock-up of the communication behavior of LULESH (3D version)
 ============================================================
 
      /  z, planes
     /
    +-----> x, cols
    |
    |
    v y, rows

  Each process/unit is responsible for a part of the overall
  simulation, called 'domain'. The size of the local domain is given
  by dx (number of columns) and dy (number of rows):




     -     E#########################F
    /     /                         /#
   /dz |---------- dx -----------| / # <-- Z1 (back plane)
  /     /                         /  #
 -  -  A#########################B <- Y0 (top plane)
    |  # --- storage order --->  #   #
    |  #                         #   H
  dy|  #                         # <----- Z0 (front plane)
    |  #                         # /
    |  #                         #/
    -  C#########################D <- Y1 (bottom plane)
       ^                         ^
       X0 (left plane)             X1 (right plane)


       A = X0Y0Z0: top-left corner (front plane)
       B = X1Y0Z0: top-right corner (front plane)
       C = X0Y1Z0: bottom-left corner (front plane)
       D = X1Y1Z0: bottom-right corner (front plane)
       E = X0Y0Z1: top-left corner (back plane)
       F = X1Y0Z1: top-right corner (back plane)
       G = X0Y1Z1: bottom-left corner (back plane)
       H = X1Y1Z1: bottom-right corner (back plane)  
       
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
  // number of local elements and nodes in each dim.
  std::array<Index_t, 3> m_nElem;
  std::array<Index_t, 3> m_nNode;

  std::vector<Real_t> m_x,  m_y,  m_z;  // position
  std::vector<Real_t> m_xd, m_yd, m_zd; // velocity
  std::vector<Real_t> m_fx, m_fy, m_fz; // force

  std::vector<Real_t> m_nodalMass;  // mass

  std::vector<Real_t> m_delv_xi;    // velocity gradient
  std::vector<Real_t> m_delv_eta;
  std::vector<Real_t> m_delv_zeta;

  Comm m_comm;
      
 public:
  Domain(std::array<int, 3> nElem,
	 std::array<int, 3> nProcs);
  
  Comm& comm() { return m_comm; }

  Index_t numProcs() const { return m_comm.numProcs(); } 
  Index_t myRank()   const { return m_comm.myRank(); }

  // number of processes in x and y direction
  Index_t px() const { return m_comm.px(); }
  Index_t py() const { return m_comm.py(); }
  Index_t pz() const { return m_comm.pz(); }

  Real_t& x(Index_t idx)  { return m_x[idx]; }
  Real_t& y(Index_t idx)  { return m_y[idx]; }
  Real_t& z(Index_t idx)  { return m_z[idx]; }
  Real_t& xd(Index_t idx) { return m_xd[idx]; }
  Real_t& yd(Index_t idx) { return m_yd[idx]; }
  Real_t& zd(Index_t idx) { return m_zd[idx]; }
  Real_t& fx(Index_t idx) { return m_fx[idx]; }
  Real_t& fy(Index_t idx) { return m_fy[idx]; }
  Real_t& fz(Index_t idx) { return m_fz[idx]; }

  Real_t& nodalMass(Index_t idx) { return m_nodalMass[idx]; }
  
  Real_t& delv_xi(Index_t idx)   { return m_delv_xi[idx]; }
  Real_t& delv_eta(Index_t idx)  { return m_delv_eta[idx]; }
  Real_t& delv_zeta(Index_t idx) { return m_delv_zeta[idx]; }

  std::array<Index_t, 3> nNode() const { return m_nNode; }
  std::array<Index_t, 3> nElem() const { return m_nElem; }
  
  Index_t nNode(int dim) const { return m_nNode[dim]; }
  Index_t nElem(int dim) const { return m_nElem[dim]; }

  Index_t maxEdgeSize() const {
    auto max = std::max(m_nNode[0], m_nNode[1]);
    max = std::max(max, m_nNode[2]);
    return max;
  }
  
  Index_t maxPlaneSize() const {
    return maxEdgeSize()*maxEdgeSize();
  }

  Real_t toValue(int c, int r, int p) {
    if( m_comm.toRank(c,r,p)<0 ) return 0.0;
    return (Real_t)p*px()*py() + (Real_t)r*px() + (Real_t)c;
  }
  
  // perform data exchange 
  void Exchange(std::array<int,3> dim,
		std::vector<Domain_member>& fields,
		Comm::Scope scope,
		Comm::Direction dir,
		Comm::Action action);

  void Print(int c, int r, int p, Domain_member f,
	     std::string str, std::array<int,3> dim);
    
  void PrintNodalMass(int col, int row, int plane);
  void PrintForce(int col, int row, int plane);
  void PrintPosVel(int col, int row, int plane);
  void PrintMonoQ(int col, int row, int plane);

  void CheckNodalMass();
};




#endif /* DOMAIN_H_INCLUDED */
