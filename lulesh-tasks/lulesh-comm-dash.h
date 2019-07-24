#ifndef LULESH_COMM_DASH_H_INCLUDED
#define LULESH_COMM_DASH_H_INCLUDED

#include <libdash.h>
#include "lulesh.h"

// forward declaration
class Domain;

#define MAX_FIELDS_PER_COMM 6

// Assume 128 byte coherence
// Assume Real_t is an "integral power of 2" bytes wide
#define CACHE_COHERENCE_PAD_REAL (128 / sizeof(Real_t))

#define CACHE_ALIGN_REAL(n)						\
  (((n) + (CACHE_COHERENCE_PAD_REAL - 1)) & ~(CACHE_COHERENCE_PAD_REAL-1))

typedef Real_t &(Domain::* Domain_member )(Index_t) ;

//
// the following is needed for the DASH version when using one-sided
// communication to find the absolute place where a plane, edge, or
// corner needs to be put. E.g., X0 is the plane where x=0, X1Z0 is
// the edge where x=1,z=0, and X1Y1Z1 is the corner at (1,1,1)
//
#define X0       0
#define X1       1
#define Y0       2
#define Y1       3
#define Z0       4
#define Z1       5
#define X0Y0     6
#define X0Y1     7
#define X1Y0     8
#define X1Y1     9
#define X0Z0    10
#define X0Z1    11
#define X1Z0    12
#define X1Z1    13
#define Y0Z0    14
#define Y0Z1    15
#define Y1Z0    16
#define Y1Z1    17
#define X0Y0Z0  18
#define X0Y0Z1  19
#define X0Y1Z0  20
#define X0Y1Z1  21
#define X1Y0Z0  22
#define X1Y0Z1  23
#define X1Y1Z0  24
#define X1Y1Z1  25

//
// additional communication helper structures that contain the stuff
// not included in the DASH version of the Domain data structure
//
class DASHComm
{
private:
  Domain& m_dom;
  dash::Array<Real_t> *m_commDataSend;
  dash::Array<Real_t> *m_commDataRecv;

public:
  typedef dash::Array<Real_t> array_type;

public:
  DASHComm(Domain& dom);
  ~DASHComm();

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
  // communication work space
  Real_t *commDataSend() { return m_commDataSend->lbegin(); }
  Real_t *commDataRecv() { return m_commDataRecv->lbegin(); }

  // given a location descriptor (e.g., X1Z0) determine its storage
  // offset
  Int_t offset(Int_t desc, Int_t xferFields);

  // determine global destination pointer for a target unit rank and
  // location descriptor
  //dash::GlobIter<Real_t, dash::Pattern<1>>
  array_type::iterator
    dest( Int_t rank, Int_t desc, Int_t xferFields );

  //dash::GlobIter<Real_t, dash::Pattern<1>>
  array_type::iterator
    src( Int_t rank, Int_t desc, Int_t xferFields );

  // 26 = 6 faces + 12 edges + 8 corners
  dash::Future<array_type::iterator> recvRequest[26];
  dash::Future<array_type::iterator> sendRequest[26];

  double wtime();
  template<typename T> T allreduce_min(T val);
  template<typename T> T reduce_max(T val);
};


void DASHCommPut(Domain& domain, DASHComm& comm,
		 Index_t xferFields, Domain_member *fieldData,
		 Index_t dx, Index_t dy, Index_t dz,
		 bool doSend, bool planeOnly);

void DASHCommSBN(Domain& domain, DASHComm& comm, int xferFields,
		 Domain_member *fieldData);

void DASHCommSyncPosVel(Domain& domain, DASHComm& comm, int xferFields,
                        Domain_member *fieldData);

void DASHCommMonoQ(Domain& domain, DASHComm& comm);


#endif /* LULESH_COMM_DASH_H_INCLUDED */
