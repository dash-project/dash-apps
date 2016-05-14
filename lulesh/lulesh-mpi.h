#ifndef LULESH_MPI_H_INCLUDED
#define LULESH_MPI_H_INCLUDED

#include <mpi.h>
#include "lulesh-dash.h"

#define MAX_FIELDS_PER_MPI_COMM 6

// Assume 128 byte coherence
// Assume Real_t is an "integral power of 2" bytes wide
#define CACHE_COHERENCE_PAD_REAL (128 / sizeof(Real_t))

#define CACHE_ALIGN_REAL(n)						\
  (((n) + (CACHE_COHERENCE_PAD_REAL - 1)) & ~(CACHE_COHERENCE_PAD_REAL-1))

typedef Real_t &(Domain::* Domain_member )(Index_t) ;

// MPI Message Tags
#define MSG_COMM_SBN      1024
#define MSG_SYNC_POS_VEL  2048
#define MSG_MONOQ         3072

// additional communication helper structures that contain the stuff
// not included in the DASH version of the Domain data structure
class Comm
{
private:
  Domain& m_dom;

public:
  Comm(Domain& dom) : m_dom(dom)
  {
    // allocate a buffer large enough for nodal ghost data
    Index_t maxPlaneSize = CACHE_ALIGN_REAL(m_dom.maxPlaneSize());
    Index_t maxEdgeSize = CACHE_ALIGN_REAL(m_dom.maxEdgeSize());

    // assume communication to 6 neighbors by default
    Index_t rowMin   = (m_dom.rowLoc() == 0)               ? 0 : 1;
    Index_t rowMax   = (m_dom.rowLoc() == m_dom.tp(1)-1)   ? 0 : 1;
    Index_t colMin   = (m_dom.colLoc() == 0)               ? 0 : 1;
    Index_t colMax   = (m_dom.colLoc() == m_dom.tp(2)-1)   ? 0 : 1;
    Index_t planeMin = (m_dom.planeLoc() == 0)             ? 0 : 1;
    Index_t planeMax = (m_dom.planeLoc() == m_dom.tp(0)-1) ? 0 : 1;

    // account for face communication
    Index_t comBufSize =
      (rowMin + rowMax + colMin + colMax + planeMin + planeMax) *
      maxPlaneSize * MAX_FIELDS_PER_MPI_COMM ;

    // account for edge communication
    comBufSize +=
      ((rowMin & colMin) + (rowMin & planeMin) + (colMin & planeMin) +
       (rowMax & colMax) + (rowMax & planeMax) + (colMax & planeMax) +
       (rowMax & colMin) + (rowMin & planeMax) + (colMin & planeMax) +
       (rowMin & colMax) + (rowMax & planeMin) + (colMax & planeMin)) *
      maxEdgeSize * MAX_FIELDS_PER_MPI_COMM ;

    // account for corner communication
    // factor of 16 is so each buffer has its own cache line
    comBufSize += ((rowMin & colMin & planeMin) +
		   (rowMin & colMin & planeMax) +
		   (rowMin & colMax & planeMin) +
		   (rowMin & colMax & planeMax) +
		   (rowMax & colMin & planeMin) +
		   (rowMax & colMin & planeMax) +
		   (rowMax & colMax & planeMin) +
		   (rowMax & colMax & planeMax)) * CACHE_COHERENCE_PAD_REAL ;

    this->commDataSend = new Real_t[comBufSize] ;
    this->commDataRecv = new Real_t[comBufSize] ;
    // prevent floating point exceptions
    memset(this->commDataSend, 0, comBufSize*sizeof(Real_t)) ;
    memset(this->commDataRecv, 0, comBufSize*sizeof(Real_t)) ;
  }

public:
  // Communication Work space
  Real_t *commDataSend;
  Real_t *commDataRecv;

  // Maximum number of block neighbors
  MPI_Request recvRequest[26]; // 6 faces + 12 edges + 8 corners
  MPI_Request sendRequest[26]; // 6 faces + 12 edges + 8 corners
};

/* doRecv flag only works with regular block structure */
void CommRecv(Domain& domain, Comm& comm, int msgType, Index_t xferFields,
              Index_t dx, Index_t dy, Index_t dz, bool doRecv, bool planeOnly);

void CommSend(Domain& domain, Comm& comm, int msgType,
              Index_t xferFields, Domain_member *fieldData,
              Index_t dx, Index_t dy, Index_t dz, bool doSend, bool planeOnly);

void CommSBN(Domain& domain, Comm& comm, int xferFields, Domain_member *fieldData);

#endif /* LULESH_MPI_H_INCLUDED */
