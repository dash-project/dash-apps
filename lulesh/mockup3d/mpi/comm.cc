
#include <iostream>
#include <assert.h>
#include <mpi.h>

#include "domain.h"
#include "comm.h"

using std::cout; using std::endl;

Comm::Comm(Domain& dom, std::array<int,3> nProcs)
  : m_dom(dom),
    m_nProcs(nProcs)
{
  MPI_Comm_rank(MPI_COMM_WORLD, &m_myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &m_numProcs);

  size_t offs=0;

  // 6 planes
  m_bufOffset[X0]=offs; offs+=MAX_FIELDS_PER_MPI_COMM*dom.maxPlaneSize();
  m_bufOffset[X1]=offs; offs+=MAX_FIELDS_PER_MPI_COMM*dom.maxPlaneSize();
  m_bufOffset[Y0]=offs; offs+=MAX_FIELDS_PER_MPI_COMM*dom.maxPlaneSize();
  m_bufOffset[Y1]=offs; offs+=MAX_FIELDS_PER_MPI_COMM*dom.maxPlaneSize();
  m_bufOffset[Z0]=offs; offs+=MAX_FIELDS_PER_MPI_COMM*dom.maxPlaneSize();
  m_bufOffset[Z1]=offs; offs+=MAX_FIELDS_PER_MPI_COMM*dom.maxPlaneSize();

  // 12 edges
  m_bufOffset[X0Y0]=offs; offs+=MAX_FIELDS_PER_MPI_COMM*dom.maxEdgeSize();
  m_bufOffset[X0Y1]=offs; offs+=MAX_FIELDS_PER_MPI_COMM*dom.maxEdgeSize();
  m_bufOffset[X1Y0]=offs; offs+=MAX_FIELDS_PER_MPI_COMM*dom.maxEdgeSize();
  m_bufOffset[X1Y1]=offs; offs+=MAX_FIELDS_PER_MPI_COMM*dom.maxEdgeSize();
  m_bufOffset[X0Z0]=offs; offs+=MAX_FIELDS_PER_MPI_COMM*dom.maxEdgeSize();
  m_bufOffset[X0Z1]=offs; offs+=MAX_FIELDS_PER_MPI_COMM*dom.maxEdgeSize();
  m_bufOffset[X1Z0]=offs; offs+=MAX_FIELDS_PER_MPI_COMM*dom.maxEdgeSize();
  m_bufOffset[X1Z1]=offs; offs+=MAX_FIELDS_PER_MPI_COMM*dom.maxEdgeSize();
  m_bufOffset[Y0Z0]=offs; offs+=MAX_FIELDS_PER_MPI_COMM*dom.maxEdgeSize();
  m_bufOffset[Y0Z1]=offs; offs+=MAX_FIELDS_PER_MPI_COMM*dom.maxEdgeSize();
  m_bufOffset[Y1Z0]=offs; offs+=MAX_FIELDS_PER_MPI_COMM*dom.maxEdgeSize();
  m_bufOffset[Y1Z1]=offs; offs+=MAX_FIELDS_PER_MPI_COMM*dom.maxEdgeSize();  

  // 8 corners
  m_bufOffset[X0Y0Z0]=offs; offs+=MAX_FIELDS_PER_MPI_COMM;
  m_bufOffset[X0Y1Z0]=offs; offs+=MAX_FIELDS_PER_MPI_COMM;
  m_bufOffset[X1Y0Z0]=offs; offs+=MAX_FIELDS_PER_MPI_COMM;
  m_bufOffset[X1Y1Z0]=offs; offs+=MAX_FIELDS_PER_MPI_COMM;
  m_bufOffset[X0Y0Z1]=offs; offs+=MAX_FIELDS_PER_MPI_COMM;
  m_bufOffset[X0Y1Z1]=offs; offs+=MAX_FIELDS_PER_MPI_COMM;
  m_bufOffset[X1Y0Z1]=offs; offs+=MAX_FIELDS_PER_MPI_COMM;
  m_bufOffset[X1Y1Z1]=offs; offs+=MAX_FIELDS_PER_MPI_COMM;  

  // offs is now total space required
  m_sendBuf.resize(offs);
  m_recvBuf.resize(offs);

  Index_t c = col();
  Index_t r = row();
  Index_t p = plane();

  // compute the rank of our neighbors in 3 dimensions
  // m_neigh[1][1][1] is the process itself (center)
  for( int x=0; x<3; ++x ) {
    for( int y=0; y<3; ++y ) {
      for( int z=0; z<3; ++z ) {
	m_neigh[x][y][z] = toRank(c-1+x, r-1+y, p-1+z);
      }
    }
  }

  // front plane
  m_trans[0][0][0] =  X0Y0Z0;
  m_trans[1][0][0] =    Y0Z0;
  m_trans[2][0][0] =  X1Y0Z0;
  m_trans[0][1][0] =    X0Z0;
  m_trans[1][1][0] =      Z0;
  m_trans[2][1][0] =    X1Z0;
  m_trans[0][2][0] =  X0Y1Z0;
  m_trans[1][2][0] =    Y1Z0;
  m_trans[2][2][0] =  X1Y1Z0;

  // middle plane
  m_trans[0][0][1] =    X0Y0;
  m_trans[1][0][1] =      Y0;
  m_trans[2][0][1] =    X1Y0;
  m_trans[0][1][1] =      X0;
  m_trans[1][1][1] =    NONE;
  m_trans[2][1][1] =      X1;
  m_trans[0][2][1] =    X0Y1;
  m_trans[1][2][1] =      Y1;
  m_trans[2][2][1] =    X1Y1;

  // back plane
  m_trans[0][0][2] =  X0Y0Z1;
  m_trans[1][0][2] =    Y0Z1;
  m_trans[2][0][2] =  X1Y0Z1;
  m_trans[0][1][2] =    X0Z1;
  m_trans[1][1][2] =      Z1;
  m_trans[2][1][2] =    X1Z1;
  m_trans[0][2][2] =  X0Y1Z1;
  m_trans[1][2][2] =    Y1Z1;
  m_trans[2][2][2] =  X1Y1Z1;
}

size_t Comm::MsgSize(std::array<int, 3> dim,
		     Transfer trans,
		     std::vector<Domain_member>& fields)
{
  size_t size=0;

  int dx = dim[0];
  int dy = dim[1];
  int dz = dim[2];
  
  switch(trans)
    {
      // planes
    case X0: case X1: size = dy*dz; break;
    case Y0: case Y1: size = dx*dz; break;
    case Z0: case Z1: size = dx*dy; break;
      
      // edges
    case X0Y0: case X0Y1: case X1Y0: case X1Y1: size = dz; break;
    case X0Z0: case X0Z1: case X1Z0: case X1Z1: size = dy; break;
    case Y0Z0: case Y0Z1: case Y1Z0: case Y1Z1: size = dx; break;
      
      // corners
    case X0Y0Z0: case X0Y0Z1: case X0Y1Z0: case X0Y1Z1:
    case X1Y0Z0: case X1Y0Z1: case X1Y1Z0: case X1Y1Z1:
      size = 1;
      break;
      
      // for completeness
    case NONE:
      size = 0;
      break;
    }
  
  return size*fields.size();
}

void Comm::Send(std::array<int, 3> dim,
		std::vector<Domain_member>& fields,
		Comm::Scope scope,
		Comm::Direction dir,
		Comm::Action action)
{
  MPI_Request reqs[26]; // 26 neighbors max
  int idx=0;

  for( int x=0; x<3; ++x ) {
    for( int y=0; y<3; ++y ) {
      for( int z=0; z<3; ++z ) {
	
	// make sure something and somewhere to send to
	if( m_neigh[x][y][z]<0 || m_trans[x][y][z]==NONE )
	  continue;
	
	// send only to higher ranks
	if( dir==Direction::One && m_neigh[x][y][z] < myRank() )
	  continue;
	
	// send planes (edges in 2D) only
	if( scope==Scope::Plane && !IS_PLANE(m_trans[x][y][z]) )
	  continue;

	PackUnpack(dim, fields, m_trans[x][y][z],
		   sendBuf(m_trans[x][y][z]), Pack::Yes, action);

	MPI_Isend(sendBuf(m_trans[x][y][z]),
		  MsgSize(dim, m_trans[x][y][z], fields),
		  MPI_DOUBLE,
		  m_neigh[x][y][z],
		  42,
		  MPI_COMM_WORLD,
		  &reqs[idx++]);
      }
    }
  }

  assert(0<= idx && idx<=26);
  
  // idx is the number of sent msgs here
  MPI_Waitall(idx, reqs, MPI_STATUSES_IGNORE);
}


void Comm::Recv(std::array<int, 3> dim,
		std::vector<Domain_member>& fields,
		Comm::Scope scope,
		Comm::Direction dir,
		Comm::Action action)
{
  for( int x=0; x<3; ++x ) {
    for( int y=0; y<3; ++y ) {
      for( int z=0; z<3; ++z ) {
	
	m_reqs[x][y][z]=MPI_REQUEST_NULL;
	
	if( m_neigh[x][y][z]<0 || m_trans[x][y][z]==NONE ) continue;
	
	// recv only from lower ranks
	if( dir==Direction::One && m_neigh[x][y][z] > myRank() )
	  continue;
	
	// recv planes (edges in 2D) only
	if( scope==Scope::Plane && !IS_PLANE(m_trans[x][y][z]) )
	  continue;

	MPI_Irecv(recvBuf(m_trans[x][y][z]),
		  MsgSize(dim, m_trans[x][y][z], fields),
		  MPI_DOUBLE,
		  m_neigh[x][y][z],
		  42,
		  MPI_COMM_WORLD,
		  &(m_reqs[x][y][z]));
      }
    }
  }
}


void Comm::Sync(std::array<int, 3> dim,
		std::vector<Domain_member>& fields,
		Comm::Scope scope,
		Comm::Direction dir,
		Comm::Action action)
{
  for( int x=0; x<3; ++x ) {
    for( int y=0; y<3; ++y ) {
      for( int z=0; z<3; ++z ) {
	
	if( m_neigh[x][y][z]<0 || m_trans[x][y][z]==NONE ) continue;
	
	// recv only from lower ranks
	if( dir==Direction::One && m_neigh[x][y][z] > myRank() )
	  continue;
	
	// recv planes (edges in 2D) only
	if( scope == Scope::Plane && !IS_PLANE(m_trans[x][y][z]) )
	  continue;
	
	MPI_Wait( &(m_reqs[x][y][z]), MPI_STATUS_IGNORE );

	PackUnpack(dim, fields, m_trans[x][y][z],
		   recvBuf(m_trans[x][y][z]), Pack::No, action);
      }
    }
  }
}


//
// pack or unpack an individual element/value
//
#define PACK_UNPACK_VAL(val_, buf_, bufidx_, pack_, action_)		\
  if( pack_==Pack::Yes ) { /* pack */					\
    buf_[bufidx_] = val_;						\
  }									\
  else { /* unpack */							\
    if( action_==Action::Add ) {					\
      val_ += buf_[bufidx_]; /* add */					\
    } else {								\
      val_ = buf_[bufidx_]; /* replace */				\
    }									\
  }

void Comm::PackUnpack(std::array<int, 3> dim,
		      std::vector<Domain_member>& fields,
		      Transfer trans,
		      Real_t *buffer,
		      Comm::Pack pack,
		      Comm::Action action)
{
  int dx = dim[0];
  int dy = dim[1];
  int dz = dim[2];

  // XY-planes - contiguous memory
  if( trans==Z0 || trans==Z1 ) {
    Index_t domidx = (trans==Z0)?0:(dx*dy*(dz-1));
    Index_t bufidx = 0;

    for( int i=0; i<dx*dy; ++i, ++domidx ) {
      assert(0<=domidx && domidx < dx*dy*dz);
      for( size_t f=0; f<fields.size(); ++f, ++bufidx ) {
	PACK_UNPACK_VAL((m_dom.*fields[f])(domidx),
			buffer, bufidx,
			pack, action);
      }
    }
  }
  
  // XZ-planes - semi-contiguous memory
  if( trans==Y0 || trans==Y1 ) {
    Index_t domidx = (trans==Y0)?0:dx*(dy-1);
    Index_t bufidx = 0;
    
    for( int i=0; i<dz; ++i, domidx+=(dx*dy)-dx ) {
      for( int j=0; j<dx; ++j, ++domidx ) {
	assert(0<=domidx && domidx < dx*dy*dz);
	for( size_t f=0; f<fields.size(); ++f, ++bufidx ) {
	  PACK_UNPACK_VAL((m_dom.*fields[f])(domidx),
			  buffer, bufidx,
			  pack, action);
	}
      }
    }
  }

  // YZ-planes - non-contiguous memory
  if( trans==X0 || trans==X1 ) {
    Index_t domidx = (trans==X0)?0:dx-1;
    Index_t bufidx = 0;
    
    for( int i=0; i<dz*dy; ++i, domidx+=dx ) {
      assert(0<=domidx && domidx < dx*dy*dz);
      for( size_t f=0; f<fields.size(); ++f, ++bufidx ) {
	PACK_UNPACK_VAL((m_dom.*fields[f])(domidx),
			buffer, bufidx,
			pack, action);
      }
    }
  }

  
  // edge in X direction
  if(trans==Y0Z0 || trans==Y0Z1 || trans==Y1Z0 || trans==Y1Z1) {
    Index_t domidx;
    switch(trans) {
    case Y0Z0:
      domidx=0; break;
    case Y1Z0:
      domidx=dx*(dy-1); break;
    case Y0Z1:
      domidx=dx*dy*(dz-1); break;
    case Y1Z1:
      domidx=dx*dy*(dz-1)+dx*(dy-1); break;
    default:
      assert(0); // should never happen
      domidx=0;
    }    
    Index_t bufidx = 0;

    // advance domidx by '1' (contiguous)
    for( int i=0; i<dx; ++i, ++domidx ) {
      assert(0<=domidx && domidx < dx*dy*dz);
      for( size_t f=0; f<fields.size(); ++f, ++bufidx ) {
	PACK_UNPACK_VAL((m_dom.*fields[f])(domidx),
			buffer, bufidx,
			pack, action);
      }
    }
  }

  // edges in Y direction
  if(trans==X0Z0 || trans==X0Z1 || trans==X1Z0 || trans==X1Z1) {
    Index_t domidx;
    switch(trans) {
    case X0Z0:
      domidx=0; break;
    case X1Z0:
      domidx=dx-1; break;
    case X0Z1:
      domidx=dx*dy*(dz-1); break;
    case X1Z1:
      domidx=dx*dy*(dz-1)+dx-1; break;
    default:
      assert(0); // should never happen
      domidx=0;
    }    
    Index_t bufidx = 0;

    // advance domidx by 'dx'
    for( int i=0; i<dy; ++i, domidx+=dx ) {
      assert( domidx < dx*dy*dz );
      for( size_t f=0; f<fields.size(); ++f, ++bufidx ) {
	PACK_UNPACK_VAL((m_dom.*fields[f])(domidx),
			buffer, bufidx,
			pack, action);
      }
    }
  }

  // edges in Z direction
  if(trans==X0Y0 || trans==X0Y1 || trans==X1Y0 || trans==X1Y1) {
    Index_t domidx;
    switch(trans) {
    case X0Y0:
      domidx=0; break;
    case X1Y0:
      domidx=dx-1; break;
    case X0Y1:
      domidx=dx*(dy-1); break;
    case X1Y1:
      domidx=dx*dy-1; break;
    default:
      assert(0); // should never happen
      domidx=0;
    }    
    Index_t bufidx = 0;

    // advance domidx by 'dx*dy' (size of XY plane)
    for( int i=0; i<dz; ++i, domidx+=dx*dy ) {
      assert( domidx < dx*dy*dz );
      for( size_t f=0; f<fields.size(); ++f, ++bufidx ) {
	PACK_UNPACK_VAL((m_dom.*fields[f])(domidx),
			buffer, bufidx,
			pack, action);
      }
    }
  }
  
  // corners
  if( IS_CORNER(trans) ) {    
    Index_t domidx;
    Index_t bufidx = 0;

    switch( trans ) {
    case X0Y0Z0: domidx = 0;            break;
    case X1Y0Z0: domidx = dx-1;         break;
    case X0Y1Z0: domidx = dx*(dy-1);    break;
    case X1Y1Z0: domidx = dx*dy-1;      break;
    case X0Y0Z1: domidx = 0         + dx*dy*(dz-1); break;
    case X1Y0Z1: domidx = dx-1      + dx*dy*(dz-1); break;
    case X0Y1Z1: domidx = dx*(dy-1) + dx*dy*(dz-1); break;
    case X1Y1Z1: domidx = dx*dy-1   + dx*dy*(dz-1); break;
    default:
      assert(0); // should never happen
      domidx=0; 
    }
    assert( domidx < dx*dy*dz );
    for( size_t f=0; f<fields.size(); ++f, ++bufidx ) {
      PACK_UNPACK_VAL((m_dom.*fields[f])(domidx),
		      buffer, bufidx,
		      pack, action);
    }
  }
}
