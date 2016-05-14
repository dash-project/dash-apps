
#include <libdash.h>
#include <iostream>
#include "lulesh-dash.h"
#include "lulesh-calc.h"

#ifdef USE_MPI
#include "lulesh-mpi.h"
#endif

using std::cout; using std::cerr; using std::endl;

Domain::Domain(const CmdLineOpts& opts) :
  // number of elements in each dim. per process
  m_nElem({{opts.nx()   ,opts.nx()   , opts.nx()  }}),

  // number of nodes in each dim. per process
  m_nNode({{opts.nx()+1 ,opts.nx()+1 , opts.nx()+1}}),

  // arrangements of processes in 3D grid
  m_ts(opts.px(), opts.py(), opts.pz()),

  // pattern for element matrix
  m_ElemPat(m_nElem[0]*opts.px(),
	    m_nElem[1]*opts.py(),
	    m_nElem[2]*opts.pz(),
	    dash::BLOCKED, dash::BLOCKED, dash::BLOCKED,
	    m_ts),

  // pattern for node matrix
  m_NodePat(m_nNode[0]*opts.px(),
	    m_nNode[1]*opts.py(),
	    m_nNode[2]*opts.pz(),
	    dash::BLOCKED, dash::BLOCKED, dash::BLOCKED,
	    m_ts)
{
  m_nodalMass.allocate(m_NodePat);
  m_x.allocate(m_NodePat);
  m_y.allocate(m_NodePat);
  m_z.allocate(m_NodePat);

  m_nodelist.allocate(m_ElemPat);
  m_volo.allocate(m_ElemPat);
  m_elemMass.allocate(m_ElemPat);

  buildMesh();

  // initialize field data
  for (Index_t i=0; i<m_ElemPat.local_size(); ++i) {
    Real_t x_local[8], y_local[8], z_local[8] ;

    Index_t *elemToNode = nodelist(i) ;
    for( Index_t lnode=0 ; lnode<8 ; ++lnode ) {
      Index_t gnode = elemToNode[lnode];

      x_local[lnode] = x(gnode);
      y_local[lnode] = y(gnode);
      z_local[lnode] = z(gnode);
    }

    // volume calculations
    Real_t volume = CalcElemVolume(x_local, y_local, z_local);

    volo(i)     = volume;
    elemMass(i) = volume;

    for (Index_t j=0; j<8; ++j) {
      Index_t idx = elemToNode[j] ;
      nodalMass(idx) += volume / Real_t(8.0) ;
    }
  }
}


Domain::~Domain()
{
  /*
  m_nodalMass.deallocate();

  m_x.deallocate();
  m_y.deallocate();
  m_z.deallocate();
  */
}


void Domain::buildMesh()
{
  auto myid = dash::myid();

  auto elem0 = m_ElemPat.extent(0);
  auto elem1 = m_ElemPat.extent(1);
  auto elem2 = m_ElemPat.extent(2);

  // get global coords for first local element
  auto gc = m_ElemPat.global({0,0,0});

  Index_t nidx = 0;
  for( Index_t plane=0; plane<m_x.local.extent(0); ++plane )
    {
      Real_t tz = Real_t(1.125)*Real_t(gc[0]+plane)/Real_t(elem0);
      for( Index_t row=0; row<m_x.local.extent(1); ++row )
	{
	  Real_t ty = Real_t(1.125)*Real_t(gc[1]+row)/Real_t(elem1);
	  for( Index_t col=0; col<m_x.local.extent(2); ++col )
	    {
	      Real_t tx = Real_t(1.125)*Real_t(gc[2]+col)/Real_t(elem2);

	      x(nidx) = tx;
	      y(nidx) = ty;
	      z(nidx) = tz;

	      ++nidx;
	    }
	}
    }

  Index_t edgeNodes = m_nodelist.local.extent(0);

  auto nplanes = m_nodelist.local.extent(0);
  auto nrows   = m_nodelist.local.extent(1);
  auto ncols   = m_nodelist.local.extent(2);

  Index_t zidx = 0; nidx=0;
  for( Index_t plane=0; plane < nplanes; ++plane ) {
    for( Index_t row=0; row < nrows; ++row ) {
      for( Index_t col=0; col < ncols; ++col )
	{
	  Index_t *localNode = nodelist(zidx);

	  localNode[0] = nidx                                      ;
	  localNode[1] = nidx                                   + 1;
	  localNode[2] = nidx                       + (ncols+1) + 1;
	  localNode[3] = nidx                       + (ncols+1)    ;
	  localNode[4] = nidx + (ncols+1)*(nrows+1)                ;
	  localNode[5] = nidx + (ncols+1)*(nrows+1)             + 1;
	  localNode[6] = nidx + (ncols+1)*(nrows+1) + (ncols+1) + 1;
	  localNode[7] = nidx + (ncols+1)*(nrows+1) + (ncols+1)    ;

	  /*
	  if( dash::myid()==0 ) {
	  cout << localNode[0] << " ";
	  cout << localNode[1] << " ";
	  cout << localNode[2] << " ";
	  cout << localNode[3] << " ";
	  cout << localNode[4] << " ";
	  cout << localNode[5] << " ";
	  cout << localNode[6] << " ";
	  cout << localNode[7] << endl;
	  }
	  */
	  
	  ++zidx;
	  ++nidx;
	}
      ++nidx;
    }
    nidx += (ncols+1);
  }
}


Real_t computeChksum(Real_t* ptr, size_t nval)
{
  Real_t sum = 0.0;

  for( size_t i=0; i<nval; ++i, ++ptr ) {
    sum += (*ptr) * Real_t(i+1)/Real_t(nval);
  }

  return sum;
}


void Domain::print_config(std::ostream& os)
{
  auto myid = dash::myid();
  os << "[ " << myid << "] colLoc()       : " << colLoc()       << std::endl;
  os << "[ " << myid << "] rowLoc()       : " << rowLoc()       << std::endl;
  os << "[ " << myid << "] planeLoc()     : " << planeLoc()     << std::endl;
  os << "[ " << myid << "] numElem()      : " << numElem()      << std::endl;
  os << "[ " << myid << "] numElem(0)     : " << numElem(0)     << std::endl;
  os << "[ " << myid << "] numElem(1)     : " << numElem(1)     << std::endl;
  os << "[ " << myid << "] numElem(2)     : " << numElem(2)     << std::endl;
  os << "[ " << myid << "] numNode()      : " << numNode()      << std::endl;
  os << "[ " << myid << "] numNode(0)     : " << numNode(0)     << std::endl;
  os << "[ " << myid << "] numNode(1)     : " << numNode(1)     << std::endl;
  os << "[ " << myid << "] numNode(2)     : " << numNode(2)     << std::endl;
  os << "[ " << myid << "] sizeX()        : " << sizeX()        << std::endl;
  os << "[ " << myid << "] sizeY()        : " << sizeY()        << std::endl;
  os << "[ " << myid << "] sizeZ()        : " << sizeZ()        << std::endl;
  os << "[ " << myid << "] maxEdgeSize()  : " << maxEdgeSize()  << std::endl;
  os << "[ " << myid << "] maxPlaneSize() : " << maxPlaneSize() << std::endl;
}
