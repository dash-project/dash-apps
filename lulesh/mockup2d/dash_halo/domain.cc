
#include <iostream>
#include <iomanip>
#include "domain.h"

Domain::Domain(const std::array<int, 2>& nElem,
	       const std::array<int, 2>& nProcs)
  : m_nElem({{nElem[0]  , nElem[1]  }}),
    m_nNode({{nElem[0]+1, nElem[1]+1}}),
    _nProcs(nProcs),
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
    stencil_spec({StencilT(-1,-1), StencilT(-1, 0), StencilT(-1, 1),
                  StencilT( 0,-1),                  StencilT( 0, 1),
                  StencilT( 1,-1), StencilT( 1, 0), StencilT( 1, 1)}),
    m_pos(m_NodePat),
    m_vel(m_NodePat),
    m_force(m_NodePat),
    m_nodalMass(m_NodePat),
    m_delv(m_ElemPat),
    halo_pos(m_pos, stencil_spec),
    halo_vel(m_vel, stencil_spec),
    halo_force(m_force, stencil_spec, CycleSpecT{Cycle::FIXED, Cycle::FIXED} ),
    halo_nodalMass(m_nodalMass, stencil_spec, CycleSpecT{Cycle::FIXED, Cycle::FIXED} ),
    halo_delv(m_delv, stencil_spec)

    //m_comm(*this, nProcs)
{
  size_t size_elem = m_nElem[0]*m_nElem[1];
  size_t size_node = m_nNode[0]*m_nNode[1];

  Index_t c = col();
  Index_t r = row();
  Real_t val = (Real_t)r * px() + (Real_t)c;

  for( size_t i=0; i<size_node; ++i ) {
    m_pos.lbegin()[i] = {val, val, val};
    m_vel.lbegin()[i] = {val, val, val};
    m_force.lbegin()[i] = {val, val, val};
    m_nodalMass.lbegin()[i] = val;
  }

  for( size_t i=0; i<size_elem; ++i ) {
    m_delv(i) = {val, val, val};
  }

  halo_nodalMass.setFixedHalos([](const std::array<dash::default_index_t,2>& coords){return 0;});

}

void Domain::exchangeNodalMass() {
  exchangeAdd(halo_nodalMass);
}

void Domain::exchangeForce() {
  exchangeAdd(halo_force);
}

template<typename HaloMatrixWrapperT>
void Domain::exchangeAdd(HaloMatrixWrapperT halo_wrapper) {
  halo_wrapper.updateHalosAsync();
  // Inner calculations
  halo_wrapper.waitHalosAsync();
  dash::barrier();
  auto it_bend = halo_wrapper.bend();
  auto& local_view = halo_wrapper.getLocalView();
  for(auto it = halo_wrapper.bbegin(); it != it_bend; ++it) {
    auto & coords = it.Coords();

    if(coords[0] == 0) {
      if(coords[1] == 0) {
        *it += it.valueAt(0) + it.valueAt(1) + it.valueAt(3);
        continue;
      }
      if(coords[1] == local_view.extent(1) - 1) {
        *it += it.valueAt(1) + it.valueAt(2) + it.valueAt(4);
        continue;
      }

      *it += it.valueAt(1);
      continue;
    }

    if(coords[0] == local_view.extent(0) - 1) {
      if(coords[1] == 0) {
        *it += it.valueAt(3) + it.valueAt(5) + it.valueAt(6);
        continue;
      }
      if(coords[1] == local_view.extent(1) - 1) {
        *it += it.valueAt(4) + it.valueAt(6) + it.valueAt(7);
        continue;
      }

      *it += it.valueAt(6);
      continue;
    }

    if(coords[1] == 0) {
      *it += it.valueAt(3);
      continue;
    }

    if(coords[1] == local_view.extent(1) - 1) {
      *it += it.valueAt(4);
      continue;
    }
  }
}

/*void Domain::Exchange(std::array<int, 2> dim,
		      std::vector<Domain_member>& fields,
		      Comm::Scope scope,
		      Comm::Direction dir,
		      Comm::Action action)
{
  // post the receives
  //m_comm.Recv(dim, fields, scope, dir, action);

  // pack data and send
  //m_comm.Send(dim, fields, scope, dir, action);

  // -----------------------------
  //   UNRELATED WORK DONE HERE
  // -----------------------------

  // wait for data transfers to finish and unpack
  //m_comm.Sync(dim, fields, scope, dir, action);
}
*/


void Domain::PrintNodalMass(int c, int r)
{
  if( col()==c && row()==r ) {
    std::cout << "------  nodalMass() of proc. "
	      << "(" << c << " " << r << ") ------ " << std::endl;
    for( int i=0; i<nNode(1); ++i ) {
      for( int j=0; j<nNode(0); ++j ) {
	std::cout << std::setw(3) << (Real_t)nodalMass(i*nNode(0)+j) << " ";
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
  if( col()==c && row()==r ) {
    std::cout << "------  fx() of proc. "
      	      << "(" << c << " " << r << ") ------ " << std::endl;
    for( int i=0; i<nNode(1); ++i ) {
      for( int j=0; j<nNode(0); ++j ) {
	std::cout << std::setw(3) << fx(i*nNode(0)+j) << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "-----------------------------------" << std::endl;
    print2dTriple(m_force, std::cout);
    std::cout << std::endl;
  }
}

/*void Domain::PrintPosVel(int c, int r)
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
*/
