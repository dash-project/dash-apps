
#include <iostream>
#include <iomanip>
#include <assert.h>
#include "domain.h"
#include "comm.h"

Domain::Domain(std::array<int, 3> nElem,
	       std::array<int, 3> nProcs)
  : m_nElem({ nElem[0]  , nElem[1]  , nElem[2]   }  ),
    m_nNode({ nElem[0]+1, nElem[1]+1, nElem[2]+1 }),
    m_comm(*this, nProcs)
{
  size_t size_elem = m_nElem[0]*m_nElem[1]*m_nElem[2];
  size_t size_node = m_nNode[0]*m_nNode[1]*m_nNode[2];

  // node-centered fields
  m_x.resize(size_node);
  m_y.resize(size_node);
  m_z.resize(size_node);
  
  m_xd.resize(size_node);
  m_yd.resize(size_node);
  m_zd.resize(size_node);

  m_fx.resize(size_node);
  m_fy.resize(size_node);
  m_fz.resize(size_node);

  m_nodalMass.resize(size_node);

  // element-centered fields
  m_delv_xi.resize(size_elem);
  m_delv_eta.resize(size_elem);
  m_delv_zeta.resize(size_elem);

  Index_t c = comm().col();
  Index_t r = comm().row();
  Index_t p = comm().plane();

  // assign value based on canonical order in
  // process grid
  Real_t val = (Real_t)p*px()*py() +
    (Real_t)r*px() + (Real_t)c;
    
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


void Domain::Exchange(std::array<int, 3> dim,
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

void Domain::Print(int c, int r, int p,
		   Domain_member f,
		   std::string str,
		   std::array<int,3> dim)
{
  if( m_comm.col()==c && m_comm.row()==r &&
      m_comm.plane()==p )
    {
      std::cout << "------  " << str << " of rank "
		<< myRank() << " "
		<< "(" << c << ", " << r << ", " << p
		<< ") ------ " << std::endl;
      
      for( int z=0; z < dim[2]; ++z ) {
	std::cout << "Plane " << z << ":" << std::endl;
	for( int y=0; y < dim[1]; ++y ) {
	  for( int x=0; x < dim[0]; ++x ) {
	    std::cout << std::setw(3)
		      << ((*this).*f)(z*dim[0]*dim[1]+
				      y*dim[0]+x)
		      << " ";
	  }
	  std::cout << std::endl;
	}
      }
      std::cout << "-----------------------------------" << std::endl;
      std::cout << std::endl;
    }
}


void Domain::PrintNodalMass(int c, int r, int p)
{
  Print(c, r, p,
	&Domain::nodalMass, "nodalMass",
	nNode() );
}


void Domain::PrintForce(int c, int r, int p)
{
  Print(c, r, p,
	&Domain::fx, "fx",
	nNode() );
}

void Domain::PrintPosVel(int c, int r, int p)
{
  Print(c, r, p,
	&Domain::x, "x",
	nNode() );
}


void Domain::PrintMonoQ(int c, int r, int p)
{
  Print(c, r, p,
	&Domain::delv_xi, "delv_xi",
	nElem() );
}

void Domain::CheckNodalMass()
{
  int c = m_comm.col();
  int r = m_comm.row();
  int p = m_comm.plane();

  int dx = m_nNode[0];
  int dy = m_nNode[1];
  int dz = m_nNode[2];

  // how many values are received from each neighbor?
  int nval[3][3][3];

  // front plane
  nval[0][0][0] =     1;
  nval[1][0][0] =    dx;
  nval[2][0][0] =     1;
  nval[0][1][0] =    dy;
  nval[1][1][0] = dx*dy;
  nval[2][1][0] =    dy;
  nval[0][2][0] =     1;
  nval[1][2][0] =    dx;
  nval[2][2][0] =     1;

  // middle plane
  nval[0][0][1] =       dz;
  nval[1][0][1] =    dx*dz;
  nval[2][0][1] =       dz;
  nval[0][1][1] =    dy*dz;
  nval[1][1][1] = dx*dy*dz;
  nval[2][1][1] =    dy*dz;
  nval[0][2][1] =       dz;
  nval[1][2][1] =    dx*dz;
  nval[2][2][1] =       dz;

  // back plane
  nval[0][0][2] =     1;
  nval[1][0][2] =    dx; 
  nval[2][0][2] =     1;
  nval[0][1][2] =    dy;
  nval[1][1][2] = dx*dy; 
  nval[2][1][2] =    dy;
  nval[0][2][2] =     1;
  nval[1][2][2] =    dx;
  nval[2][2][2] =     1;
  
  
  // determine actual sum over whole domain
  Real_t actual=0.0;
  for( int i=0; i<dx*dy*dz; ++i ) {
    actual += nodalMass(i);
  }

  // compute expected sum over domain by accounting
  // for all contributions from neighboring processes
  Real_t expected=0.0;
  
  for( int i=-1; i<2; ++i ) {
    for( int j=-1; j<2; ++j ) {
      for( int k=-1; k<2; ++k ) {
	if( m_comm.toRank(c+i,r+j,p+k)<0 )
	  continue;

	// neighbor's assigned value
	Real_t val = toValue(c+i,r+j,p+k);
	expected += nval[1+i][1+j][1+k] * val;
      }
    }
  }
  
  if( expected!=actual ) {
    std::cout << std::flush;
    std::cout << "CheckNodalMass failed on Rank "
	      << myRank()
	      << " expected: " << expected
	      << " actual: " << actual << std::endl;
    std::cout << std::flush;
  }
}
