
// for allreduce -- FIXME
#include <mpi.h>

#include "lulesh.h"
#include "lulesh-dash.h"
#include "lulesh-comm-dash.h"

DASHComm::DASHComm(Domain& dom) : m_dom(dom)
{
  // allocate a buffer large enough for nodal ghost data
  Index_t maxPlaneSize = CACHE_ALIGN_REAL( dom.maxPlaneSize() );
  Index_t maxEdgeSize  = CACHE_ALIGN_REAL( dom.maxEdgeSize()  );

  // account for face communication (6 faces)
  Index_t comBufSize = 6 * maxPlaneSize * MAX_FIELDS_PER_COMM;

  // account for edge communication (12 edges)
  comBufSize += 12 * maxEdgeSize * MAX_FIELDS_PER_COMM;

  // account for corner communication (8 corners)
  // factor of 16 is so each buffer has its own cache line
  comBufSize +=  8 * CACHE_COHERENCE_PAD_REAL;

  m_commDataSend = new dash::Array<Real_t>(dash::size()*comBufSize,
					   dash::BLOCKED);

  m_commDataRecv = new dash::Array<Real_t>(dash::size()*comBufSize,
					   dash::BLOCKED);

  memset( m_commDataSend->lbegin(), -1, comBufSize*sizeof(Real_t) );
  memset( m_commDataRecv->lbegin(), -1, comBufSize*sizeof(Real_t) );

  m_commDataRecv->barrier();
}

DASHComm::~DASHComm()
{
  delete m_commDataSend;
  delete m_commDataRecv;
}

void DASHComm::ExchangeNodalMass()
{
  Domain_member fieldData;
  fieldData = &Domain::nodalMass;

  DASHComm& comm = *this;
  Domain&  dom  = m_dom;

  // Initial domain boundary communication
  DASHCommPut(dom, comm, 1, &fieldData,
	  dom.sizeX() + 1, dom.sizeY() + 1, dom.sizeZ() +  1,
	  true, false);

  dash::tasks::async_barrier();

  //dash::barrier();
  DASHCommSBN(dom, comm, 1, &fieldData);

  dash::tasks::complete();
}

void DASHComm::Recv_PosVel()
{
  // empty
}

void DASHComm::Send_PosVel()
{
  //Domain_member fieldData[6];

  static Domain_member *fieldData = nullptr;
  if (fieldData == nullptr) {
    fieldData    = new Domain_member[6];
    fieldData[0] = &Domain::x;
    fieldData[1] = &Domain::y;
    fieldData[2] = &Domain::z;
    fieldData[3] = &Domain::xd;
    fieldData[4] = &Domain::yd;
    fieldData[5] = &Domain::zd;
  }

  DASHComm& comm = *this;
  Domain&  dom  = m_dom;

  DASHCommPut(dom, comm, 6, fieldData,
	  dom.sizeX() + 1, dom.sizeY() + 1, dom.sizeZ() + 1,
	  false, false);
}

void DASHComm::Sync_PosVel()
{
  DASHComm& comm = *this;
  Domain&  dom  = m_dom;
  //Domain_member fieldData[6];
  static Domain_member *fieldData = nullptr;
  if (fieldData == nullptr) {
    fieldData    = new Domain_member[6];
    fieldData[0] = &Domain::x;
    fieldData[1] = &Domain::y;
    fieldData[2] = &Domain::z;
    fieldData[3] = &Domain::xd;
    fieldData[4] = &Domain::yd;
    fieldData[5] = &Domain::zd;
  }
//
  dash::tasks::async_barrier();

  // dummy task to emulate CONCURRENT dependency
  dash::tasks::async(
    "Sync_PosVel_DUMMY",
    [](){
      // nothing to do
    },
    dash::tasks::out(dom.x(0))
  );


  DASHCommSyncPosVel(dom, comm, 6, fieldData);

  // dummy task to emulate CONCURRENT dependency
  dash::tasks::async(
    "Sync_PosVel_DUMMY2",
    [](){
      // nothing to do
    },
    dash::tasks::out(dom.x(0))
  );
  dash::tasks::async_barrier();
}


void DASHComm::Recv_Force()
{
  // empty
}

void DASHComm::Send_Force()
{
  DASHComm& comm = *this;
  Domain&  dom  = m_dom;
  //Domain_member fieldData[3];
  static Domain_member *fieldData = nullptr;
  if (fieldData == nullptr) {
    fieldData    = new Domain_member[3];
    fieldData[0] = &Domain::fx;
    fieldData[1] = &Domain::fy;
    fieldData[2] = &Domain::fz;
  }

  DASHCommPut(dom, comm, 3, fieldData,
              dom.sizeX() + 1, dom.sizeY() + 1, dom.sizeZ() +  1,
              true, false);
}

void DASHComm::Sync_Force()
{
  //Domain_member fieldData[3];
  static Domain_member *fieldData = nullptr;
  if (fieldData == nullptr) {
    fieldData    = new Domain_member[3];
    fieldData[0] = &Domain::fx;
    fieldData[1] = &Domain::fy;
    fieldData[2] = &Domain::fz;
  }

  DASHComm& comm = *this;
  Domain&  dom  = m_dom;

  dash::tasks::async_barrier();

  // dummy task to emulate CONCURRENT dependency
  dash::tasks::async(
    "Sync_Force_DUMMY",
    [](){
      // nothing to do
    },
    dash::tasks::out(dom.fx(0))
  );

  DASHCommSBN(dom, comm, 3, fieldData);

  // dummy task to emulate CONCURRENT dependency
  dash::tasks::async(
    "Sync_Force_DUMMY2",
    [](){
      // nothing to do
    },
    dash::tasks::out(dom.fx(0))
  );
  dash::tasks::async_barrier();
}

void DASHComm::Recv_MonoQ()
{
  // empty
}

void DASHComm::Send_MonoQ()
{
  DASHComm& comm = *this;
  Domain&  dom  = m_dom;
  //Domain_member fieldData[3];
  // Transfer veloctiy gradients in the first order elements
  // problem->commElements->Transfer(CommElements::monoQ);
  static Domain_member *fieldData = nullptr;
  if (fieldData == nullptr) {
    fieldData    = new Domain_member[3];
    fieldData[0] = &Domain::delv_xi;
    fieldData[1] = &Domain::delv_eta;
    fieldData[2] = &Domain::delv_zeta;
  }

  DASHCommPut(dom, comm, 3, fieldData,
	  dom.sizeX(), dom.sizeY(), dom.sizeZ(),
	  true, true);
}

void DASHComm::Sync_MonoQ()
{
  DASHComm& comm = *this;
  Domain&  dom  = m_dom;

  dash::tasks::async_barrier();

  // dummy task to emulate CONCURRENT dependency
  dash::tasks::async(
    "Sync_MonoQ_DUMMY",
    [](){
      // nothing to do
    },
    dash::tasks::out(dom.delv_xi(0))
  );

  DASHCommMonoQ(dom,comm);

  // dummy task to emulate CONCURRENT dependency
  dash::tasks::async(
    "Sync_MonoQ_DUMMY2",
    [](){
      // nothing to do
    },
    dash::tasks::out(dom.delv_xi(0))
  );
}

Int_t DASHComm::offset(Int_t desc, Int_t xferFields)
{
  Index_t maxPlaneSize = CACHE_ALIGN_REAL( m_dom.maxPlaneSize() );
  Index_t maxEdgeSize  = CACHE_ALIGN_REAL( m_dom.maxEdgeSize()  );

  Int_t pmsg = desc >=  6 ?  6 : desc; desc -= pmsg;
  Int_t emsg = desc >= 12 ? 12 : desc; desc -= emsg;
  Int_t cmsg = desc;

  assert(0 <= pmsg && pmsg <= 6);
  assert(0 <= emsg && emsg <= 12);
  assert(0 <= cmsg && cmsg <= 8);

  auto offs =
    (pmsg * maxPlaneSize +
     emsg * maxEdgeSize) * MAX_FIELDS_PER_COMM +
    cmsg * CACHE_COHERENCE_PAD_REAL;

  return offs;
}

dash::GlobIter<Real_t, dash::Pattern<1>>
DASHComm::dest(Int_t rank, Int_t desc, Int_t xferFields)
{
  auto& pat = m_commDataRecv->pattern();

  auto offs = offset(desc, xferFields);

  dash::GlobIter<Real_t, dash::Pattern<1> > it
    = m_commDataRecv->begin() + pat.global_index(dash::team_unit_t(rank),
                                                 {offs});

  return it;
}

dash::GlobIter<Real_t, dash::Pattern<1>>
DASHComm::src(Int_t rank, Int_t desc, Int_t xferFields)
{
  auto& pat = m_commDataSend->pattern();

  auto offs = offset(desc, xferFields);

  dash::GlobIter<Real_t, dash::Pattern<1> > it
    = m_commDataSend->begin() + pat.global_index(dash::team_unit_t(rank),
                                                 {offs});

  return it;
}


double DASHComm::wtime()
{
  // FIXME
  return MPI_Wtime();
}

template<>
double DASHComm::allreduce_min<double>(double val)
{
  double res;
  // FIXME
  MPI_Request request;
  MPI_Iallreduce(&val, &res, 1, MPI_DOUBLE, MPI_MIN,
		MPI_COMM_WORLD, &request);
  do {
    int flag;
    MPI_Test(&request, &flag, MPI_STATUS_IGNORE);
    if (flag) break;
    dash::tasks::yield();
  } while (1);
  return res;
}

template<>
double DASHComm::reduce_max<double>(double val)
{
  double res;
  // FIXME
  MPI_Reduce(&val, &res, 1, MPI_DOUBLE, MPI_MAX, 0,
	     MPI_COMM_WORLD);
  return res;
}


