
#include <iostream>
#include <mpi.h>

#include "domain.h"
#include "comm.h"

using std::cout; using std::endl;

Comm::Comm(Domain& dom, std::array<int,2> nProcs)
  : m_dom(dom),
    m_nProcs(nProcs),
    m_trans{{ X0Y0, Y0  , X1Y0 },
            { X0  , NONE, X1   },
            { X0Y1, Y1  , X1Y1 }},
    m_neigh{{-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}}
{
  size_t offs=0;

  MPI_Comm_rank(MPI_COMM_WORLD, &m_myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &m_numProcs);
  
  // 4 edges
  m_bufOffset[X0]=offs; offs+=MAX_FIELDS_PER_MPI_COMM*dom.maxEdgeSize();
  m_bufOffset[X1]=offs; offs+=MAX_FIELDS_PER_MPI_COMM*dom.maxEdgeSize();
  m_bufOffset[Y0]=offs; offs+=MAX_FIELDS_PER_MPI_COMM*dom.maxEdgeSize();
  m_bufOffset[Y1]=offs; offs+=MAX_FIELDS_PER_MPI_COMM*dom.maxEdgeSize();

  // 4 corners
  m_bufOffset[X0Y0]=offs; offs+=MAX_FIELDS_PER_MPI_COMM;
  m_bufOffset[X0Y1]=offs; offs+=MAX_FIELDS_PER_MPI_COMM;
  m_bufOffset[X1Y0]=offs; offs+=MAX_FIELDS_PER_MPI_COMM;
  m_bufOffset[X1Y1]=offs; offs+=MAX_FIELDS_PER_MPI_COMM;

  // offs is now total space required
  m_sendBuf.resize(offs);
  m_recvBuf.resize(offs);

  Index_t c = col();
  Index_t r = row();
  
  m_neigh[0][0] = toRank(c-1, r-1);
  m_neigh[0][1] = toRank(c,   r-1);
  m_neigh[0][2] = toRank(c+1, r-1);

  m_neigh[1][0] = toRank(c-1, r);
  m_neigh[1][1] = toRank(c,   r);
  m_neigh[1][2] = toRank(c+1, r);

  m_neigh[2][0] = toRank(c-1, r+1);
  m_neigh[2][1] = toRank(c,   r+1);
  m_neigh[2][2] = toRank(c+1, r+1);
}

size_t Comm::MsgSize(std::array<int, 2> dim,
		     Transfer trans,
		     std::vector<Domain_member>& fields)
{
  size_t size=0;

  int dx = dim[0];
  int dy = dim[1];
  
  switch(trans) {
    // edges in y direction
  case X0: case X1:
    size = dy*fields.size();
    break;
    
    // edges in x direction
  case Y0: case Y1:
    size = dx*fields.size();
    break;
    
    // corners
  case X0Y0: case X0Y1: case X1Y0: case X1Y1:
    size = fields.size();
    break;

    // for completeness
  case NONE:
    size = 0;
    break;
  }
  return size;
}

void Comm::Send(std::array<int, 2> dim,
		std::vector<Domain_member>& fields,
		Comm::Scope scope,
		Comm::Direction dir,
		Comm::Action action)
{
  MPI_Request reqs[8]; // 8 neighbors max
  int idx=0;
  
  for( int x=0; x<3; ++x ) {
    for( int y=0; y<3; ++y ) {

      // make sure something and somewhere to send to
      if( m_neigh[x][y]<0 || m_trans[x][y]==NONE )
	continue;

      // send only to higher ranks
      if( dir==Direction::One && m_neigh[x][y] < myRank() )
	continue;

      // send planes (edges in 2D) only
      if( scope==Scope::Plane && !(m_trans[x][y]==X0 ||
				   m_trans[x][y]==X1 ||
				   m_trans[x][y]==Y0 ||
				   m_trans[x][y]==Y1) )
	continue;
			
      PackUnpack(dim, fields, m_trans[x][y],
		 sendBuf(m_trans[x][y]), Pack::Yes, action);

      MPI_Isend(sendBuf(m_trans[x][y]),
		MsgSize(dim, m_trans[x][y], fields),
		MPI_DOUBLE,
		m_neigh[x][y],
		42,
		MPI_COMM_WORLD,
		&reqs[idx++]);
    }
  }

  // idx is the number of sent msgs here
  MPI_Waitall(idx, reqs, MPI_STATUSES_IGNORE);
}


void Comm::Recv(std::array<int, 2> dim,
		std::vector<Domain_member>& fields,
		Comm::Scope scope,
		Comm::Direction dir,
		Comm::Action action)
{
  for( int x=0; x<3; ++x ) {
    for( int y=0; y<3; ++y ) {
      m_reqs[x][y]=MPI_REQUEST_NULL;
	
      if( m_neigh[x][y]<0 || m_trans[x][y]==NONE ) continue;

      // recv only from lower ranks
      if( dir==Direction::One && m_neigh[x][y] > myRank() )
	continue;

      // recv planes (edges in 2D) only
      if( scope==Scope::Plane && !(m_trans[x][y]==X0 ||
				   m_trans[x][y]==X1 ||
				   m_trans[x][y]==Y0 ||
				   m_trans[x][y]==Y1) )
	continue;

      MPI_Irecv(recvBuf(m_trans[x][y]),
		MsgSize(dim, m_trans[x][y], fields),
		MPI_DOUBLE,
		m_neigh[x][y],
		42,
		MPI_COMM_WORLD,
		&(m_reqs[x][y]));
    }
  }
}


void Comm::Sync(std::array<int, 2> dim,
		std::vector<Domain_member>& fields,
		Comm::Scope scope,
		Comm::Direction dir,
		Comm::Action action)
{
  for( int x=0; x<3; ++x ) {
    for( int y=0; y<3; ++y ) {
      if( m_neigh[x][y]<0 || m_trans[x][y]==NONE ) continue;

      // recv only from lower ranks
      if( dir==Direction::One && m_neigh[x][y] > myRank() )
	continue;

      // recv planes (edges in 2D) only
      if( scope == Scope::Plane && !(m_trans[x][y]==X0 ||
				     m_trans[x][y]==X1 ||
				     m_trans[x][y]==Y0 ||
				     m_trans[x][y]==Y1) )
	continue;
      
      MPI_Wait( &(m_reqs[x][y]), MPI_STATUS_IGNORE );

      PackUnpack(dim, fields, m_trans[x][y],
		 recvBuf(m_trans[x][y]), Pack::No, action);
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

void Comm::PackUnpack(std::array<int, 2> dim,
		      std::vector<Domain_member>& fields,
		      Transfer trans,
		      Real_t *buffer,
		      Comm::Pack pack,
		      Comm::Action action)
{
  int dx = dim[0];
  int dy = dim[1];

  // edge - contiguous memory
  if( trans==Y0 || trans==Y1 ) {
    Index_t domidx = (trans==Y0)?0:(dx*(dy-1));
    Index_t bufidx = 0;

    for( int i=0; i<dx; ++i, ++domidx ) {
      for( size_t f=0; f<fields.size(); ++f, ++bufidx ) {
	PACK_UNPACK_VAL((m_dom.*fields[f])(domidx),
			buffer, bufidx,
			pack, action);
      }
    }
  }

  // edge - non-contiguous memory
  if( trans==X0 || trans==X1 ) {
    Index_t domidx = (trans==X0)?0:dx-1;
    Index_t bufidx = 0;

    for( int i=0; i<dy; ++i, domidx+=dx ) {
      for( size_t f=0; f<fields.size(); ++f, ++bufidx ) {
	PACK_UNPACK_VAL((m_dom.*fields[f])(domidx),
			buffer, bufidx,
			pack, action);
      }
    }
  }

  // corners
  if( trans==X0Y0 || trans==X0Y1 ||
      trans==X1Y0 || trans==X1Y1 ) {

    Index_t domidx;
    Index_t bufidx = 0;
    
    if( trans==X0Y0 ) domidx = 0;
    if( trans==X1Y0 ) domidx = dx-1;
    if( trans==X0Y1 ) domidx = dx*(dy-1);
    if( trans==X1Y1 ) domidx = dx*dy-1;
    
    for( size_t f=0; f<fields.size(); ++f, ++bufidx ) {
      PACK_UNPACK_VAL((m_dom.*fields[f])(domidx),
		      buffer, bufidx,
		      pack, action);
    }
  }
}
