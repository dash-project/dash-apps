
#include <iostream>
#include <iomanip>
#include "domain.h"
#include "comm.h"

Domain::Domain(std::array<int, 2> nElem,
	       std::array<int, 2> nProcs)
  : m_nElem({{nElem[0]  , nElem[1]  }}),
    m_nNode({{nElem[0]+1, nElem[1]+1}}),
    
    // process arrangement
    m_teamspec(nProcs[0], nProcs[1]),
    
    // pattern for element matrix
    m_ElemPat(nElem[0]*nProcs[0],
	      nElem[1]*nProcs[1],
	      dash::BLOCKED, dash::BLOCKED,
	      m_teamspec),
    
    // pattern for node matrix
    m_NodePat((nElem[0]+1)*nProcs[0],
	      (nElem[1]+1)*nProcs[1],
	      dash::BLOCKED, dash::BLOCKED,
	      m_teamspec),
    
    m_comm(*this, nProcs)
{
  size_t size_elem = m_nElem[0]*m_nElem[1];
  size_t size_node = m_nNode[0]*m_nNode[1];

  m_x.allocate(m_NodePat);
  m_y.allocate(m_NodePat);
  m_z.allocate(m_NodePat);
  
  m_xd.allocate(m_NodePat);
  m_yd.allocate(m_NodePat);
  m_zd.allocate(m_NodePat);

  m_fx.allocate(m_NodePat);
  m_fy.allocate(m_NodePat);
  m_fz.allocate(m_NodePat);

  m_nodalMass.allocate(m_NodePat);

  m_delv_xi.allocate(m_ElemPat);
  m_delv_eta.allocate(m_ElemPat);
  m_delv_zeta.allocate(m_ElemPat);

  Index_t c = comm().col();
  Index_t r = comm().row();
  Real_t val = (Real_t)r * px() + (Real_t)c;
    
  for( size_t i=0; i<size_node; ++i ) {
    x(i)  = y(i)  = z(i)  = val;
    xd(i) = yd(i) = zd(i) = val;
    fx(i) = fy(i) = fz(i) = val;
    nodalMass(i) = val;
  }

  for( size_t i=0; i<size_elem; ++i ) {
    delv_xi(i)   = val;
    delv_eta(i)  = val;
    delv_zeta(i) = val;
  }
}


void Domain::Exchange(std::array<int, 2> dim,
		      std::vector<Domain_member>& fields,
		      Comm::Scope scope,
		      Comm::Direction dir,
		      Comm::Action action)
{
  // post the receives
  m_comm.Recv(dim, fields, scope, dir, action);

  // pack data and send
  m_comm.Send(dim, fields, scope, dir, action);

  // -----------------------------
  //   UNRELATED WORK DONE HERE
  // -----------------------------
  
  // wait for data transfers to finish and unpack
  m_comm.Sync(dim, fields, scope, dir, action);
}


void Domain::PrintNodalMass(int c, int r)
{
  if( m_comm.col()==c && m_comm.row()==r ) {
    std::cout << "------  nodalMass() of proc. "
	      << "(" << c << " " << r << ") ------ " << std::endl;
    for( int i=0; i<nNode(1); ++i ) {
      for( int j=0; j<nNode(0); ++j ) {
	std::cout << std::setw(3) << nodalMass(i*nNode(0)+j) << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "-----------------------------------" << std::endl;
    print2d(m_nodalMass, std::cout);
    std::cout << std::endl;
  }
}


void Domain::PrintForce(int c, int r)
{
  if( m_comm.col()==c && m_comm.row()==r ) {
    std::cout << "------  fx() of proc. "
      	      << "(" << c << " " << r << ") ------ " << std::endl;
    for( int i=0; i<nNode(1); ++i ) {
      for( int j=0; j<nNode(0); ++j ) {
	std::cout << std::setw(3) << fx(i*nNode(0)+j) << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "-----------------------------------" << std::endl;
    print2d(m_fx, std::cout);
    std::cout << std::endl;
  }
}

void Domain::PrintPosVel(int c, int r)
{
  if( m_comm.col()==c && m_comm.row()==r ) {
    std::cout << "------  x() of proc. "
	      << "(" << c << " " << r << ") ------ " << std::endl;
      for( int i=0; i<nNode(1); ++i ) {
      for( int j=0; j<nNode(0); ++j ) {
	std::cout << std::setw(3) << x(i*nNode(0)+j) << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "-----------------------------------" << std::endl;
    print2d(m_x, std::cout);
    std::cout << std::endl;
  }
}


void Domain::PrintMonoQ(int c, int r)
{
  if( m_comm.col()==c && m_comm.row()==r ) {
    std::cout << "------  delv_xi() of proc. "
	      << "(" << c << " " << r << ") ------ " << std::endl;
    for( int i=0; i<nElem(1); ++i ) {
      for( int j=0; j<nElem(0); ++j ) {
	std::cout << std::setw(3) << delv_xi(i*nElem(0)+j) << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "-----------------------------------" << std::endl;
    print2d(m_delv_xi, std::cout);
    std::cout << std::endl;
  }
}
