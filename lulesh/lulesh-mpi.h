#ifndef LULESH_MPI_H_INCLUDED
#define LULESH_MPI_H_INCLUDED

#include "lulesh.h"
#include <mpi.h>

// forward declaration
class Domain;

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

//
// additional communication helper structures that contain the stuff
// not included in the DASH version of the Domain data structure
//
class MPIComm
{
private:
  Domain& m_dom;

public:
  MPIComm(Domain& dom);
  ~MPIComm();

  void ExchangeNodalMass(); // 1 field: nodalMass

  void Recv_PosVel();       // 6 fields: (x,y,z), (dx,dy,dz)
  void Send_PosVel();
  void Sync_PosVel();

  void Recv_Force();        // 3 fields: (fx,fy,fz)
  void Send_Force();
  void Sync_Force();

  void Recv_MonoQ();        // 3 fields: delv_xi, delv_eta, delv_zeta
  void Send_MonoQ();
  void Sync_MonoQ();

public:
  // Communication Work space
  Real_t *commDataSend;
  Real_t *commDataRecv;

  // Maximum number of block neighbors
  MPI_Request recvRequest[26]; // 6 faces + 12 edges + 8 corners
  MPI_Request sendRequest[26]; // 6 faces + 12 edges + 8 corners
};


// doRecv flag only works with regular block structure
void CommRecv(Domain& domain, MPIComm& comm, int msgType,
	      Index_t xferFields, Index_t dx, Index_t dy, Index_t dz,
	      bool doRecv, bool planeOnly);

void CommSend(Domain& domain, MPIComm& comm, int msgType,
              Index_t xferFields, Domain_member *fieldData,
              Index_t dx, Index_t dy, Index_t dz,
	      bool doSend, bool planeOnly);

void CommSBN(Domain& domain, MPIComm& comm, int xferFields,
	     Domain_member *fieldData);

void CommSyncPosVel(Domain& domain, MPIComm& comm);

void CommMonoQ(Domain& domain, MPIComm& comm);


#endif /* LULESH_MPI_H_INCLUDED */
