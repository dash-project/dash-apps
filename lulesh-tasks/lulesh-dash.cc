
#include <libdash.h>
#include <iostream>

#include "lulesh-dash.h"
#include "lulesh-calc.h"
#include "lulesh-util.h"

#ifdef USE_MPI
#include <mpi.h>
#include "lulesh-comm-mpi.h"
#endif

#ifdef USE_DASH
#include "lulesh-comm-dash.h"
#endif

#include "extrae.h"

#define DASH_TASKLOOP_FACTOR 6

using std::cout; using std::cerr; using std::endl;

template <typename T>
static inline
T *Allocate(size_t size)
{
  return static_cast<T *>(malloc(sizeof(T)*size)) ;
}

template <typename T>
static inline
void Release(T **ptr)
{
  if (*ptr != NULL) {
    free(*ptr) ;
    *ptr = NULL ;
  }
}

template<>
double Domain::allreduce_min<double>(double val)
{
  return m_comm.allreduce_min(val);
}

template<>
double Domain::reduce_max<double>(double val)
{
  return m_comm.reduce_max(val);
}


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
	    m_ts),

  m_region(opts.numReg(), opts.cost(),
	   opts.balance(), numElem() ),

  m_comm(*this)
{
  //
  // node - centered fields
  //
#if 0
  m_x.allocate(m_NodePat);
  m_y.allocate(m_NodePat);
  m_z.allocate(m_NodePat);

  m_xd.allocate(m_NodePat);
  m_yd.allocate(m_NodePat);
  m_zd.allocate(m_NodePat);

  m_xdd.allocate(m_NodePat);
  m_ydd.allocate(m_NodePat);
  m_zdd.allocate(m_NodePat);

  m_fx.allocate(m_NodePat);
  m_fy.allocate(m_NodePat);
  m_fz.allocate(m_NodePat);

  m_nodalMass.allocate(m_NodePat);
#endif

  m_x.resize(numNode());
  m_y.resize(numNode());
  m_z.resize(numNode());

  m_xd.resize(numNode());
  m_yd.resize(numNode());
  m_zd.resize(numNode());

  m_xdd.resize(numNode());
  m_ydd.resize(numNode());
  m_zdd.resize(numNode());

  m_fx.resize(numNode());
  m_fy.resize(numNode());
  m_fz.resize(numNode());

  m_nodalMass.resize(numNode());

  //
  // element-centered fields
  //
  m_e.resize(numElem());
  m_p.resize(numElem());
  m_q.resize(numElem());
  m_ql.resize(numElem());
  m_qq.resize(numElem());

  m_v.resize(numElem());
  m_volo.resize(numElem());
  m_delv.resize(numElem());
  m_vdov.resize(numElem());
  m_elemMass.resize(numElem());

  m_nodelist.resize(numElem());

  AllocateStrains(numElem());

  Int_t allElem = numElem() +           /* local elem */
    2*sizeX()*sizeY() + /* plane ghosts */
    2*sizeX()*sizeZ() + /* row ghosts */
    2*sizeY()*sizeZ() ; /* col ghosts */

  AllocateGradients(numElem(), allElem);

  m_arealg.resize(numElem());
  m_ss.resize(numElem());

  m_vnew.resize(numElem());

  Index_t edgeElems = opts.nx();
  Index_t edgeNodes = edgeElems+1;

  BuildMesh();
  SetupThreadSupportStructures();

  // Setup region index sets. For now, these are constant sized
  // throughout the run, but could be changed every cycle to
  // simulate effects of ALE on the lagrange solver


  // Setup symmetry nodesets
  SetupSymmetryPlanes(edgeNodes);

  // Setup element connectivities
  SetupElementConnectivities(edgeElems);

  // Setup symmetry planes and free surface boundary arrays
  SetupBoundaryConditions(edgeElems);

  // Initialize field data
  InitializeFieldData();

  // Deposit initial energy
  DepositInitialEnergy(opts.nx());
}


Domain::~Domain()
{
  /* XXX
  m_nodalMass.deallocate();

  m_x.deallocate();
  m_y.deallocate();
  m_z.deallocate();
  */
}


void Domain::BuildMesh()
{
  auto elem0 = m_ElemPat.extent(0);
  auto elem1 = m_ElemPat.extent(1);
  auto elem2 = m_ElemPat.extent(2);

  // get global coords for first local element
  auto gc = m_ElemPat.global({0,0,0});

  Index_t nidx = 0;
  for( Index_t plane=0; plane<m_NodePat.local_extent(0); ++plane )
  {
    Real_t tz = Real_t(1.125)*Real_t(gc[0]+plane)/Real_t(elem0);
    for( Index_t row=0; row<m_NodePat.local_extent(1); ++row )
    {
      Real_t ty = Real_t(1.125)*Real_t(gc[1]+row)/Real_t(elem1);
      for( Index_t col=0; col<m_NodePat.local_extent(2); ++col )
      {
        Real_t tx = Real_t(1.125)*Real_t(gc[2]+col)/Real_t(elem2);

        x(nidx) = tx;
        y(nidx) = ty;
        z(nidx) = tz;

        ++nidx;
      }
    }
  }

  auto nplanes = m_ElemPat.local_extent(0);
  auto nrows   = m_ElemPat.local_extent(1);
  auto ncols   = m_ElemPat.local_extent(2);

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

void Domain::InitialBoundaryExchange()
{
  m_comm.ExchangeNodalMass();
}


void Domain::TimeIncrement()
{
  Domain& domain = (*this);
  dash::tasks::ASYNC(
    [=, &domain](){
      EXTRAE_ENTER(TIME_INCREMENT);
      Real_t targetdt = domain.stoptime() - domain.time() ;

      if ((domain.dtfixed() <= Real_t(0.0)) && (domain.cycle() != Int_t(0))) {
        Real_t ratio ;
        Real_t olddt = domain.deltatime() ;

        /* This will require a reduction in parallel */
        Real_t gnewdt = Real_t(1.0e+20) ;
        Real_t newdt ;
        if (domain.dtcourant() < gnewdt) {
          gnewdt = domain.dtcourant() / Real_t(2.0) ;
        }
        if (domain.dthydro() < gnewdt) {
          gnewdt = domain.dthydro() * Real_t(2.0) / Real_t(3.0) ;
        }

        newdt = m_comm.allreduce_min(gnewdt);

        ratio = newdt / olddt ;
        if (ratio >= Real_t(1.0)) {
          if (ratio < domain.deltatimemultlb()) {
            newdt = olddt ;
          }
          else if (ratio > domain.deltatimemultub()) {
            newdt = olddt*domain.deltatimemultub() ;
          }
        }

        if (newdt > domain.dtmax()) {
          newdt = domain.dtmax() ;
        }
        domain.deltatime() = newdt ;
      }

      /* TRY TO PREVENT VERY SMALL SCALING ON THE NEXT CYCLE */
      if ((targetdt > domain.deltatime()) &&
          (targetdt < (Real_t(4.0) * domain.deltatime() / Real_t(3.0))) ) {
        targetdt = Real_t(2.0) * domain.deltatime() / Real_t(3.0) ;
      }

      if (targetdt < domain.deltatime()) {
        domain.deltatime() = targetdt ;
      }

      domain.time() += domain.deltatime() ;

      ++domain.cycle() ;
      EXTRAE_EXIT(TIME_INCREMENT);
    },
    DART_PRIO_HIGH,
    dash::tasks::in( domain.dtcourant()),
    dash::tasks::in( domain.dthydro()),
    dash::tasks::out(domain.deltatime())
  );
}


void Domain::LagrangeLeapFrog()
{
  LagrangeNodal();
  LagrangeElements();

  CalcTimeConstraintsForElems(*this);

}

void Domain::LagrangeNodal()
{
  Real_t ucut = u_cut();

  CalcForceForNodes();

  dash::tasks::ASYNC(
    [=](){
      EXTRAE_ENTER(CALCACCELERATIONFORNODES);
      CalcAccelerationForNodes(numNode());
      ApplyAccelerationBoundaryConditionsForNodes();
      EXTRAE_EXIT(CALCACCELERATIONFORNODES);
    },
    dash::tasks::in( this->fx(0)),
    dash::tasks::out(this->xdd(0))
  );
  dash::tasks::ASYNC(
    [=](){
      EXTRAE_ENTER(CALCVELPOSFORNODES);
      const Real_t delt = deltatime();
      CalcVelocityForNodes(delt, ucut, numNode()) ;
      CalcPositionForNodes(delt, numNode());
      dash::tasks::complete();
      EXTRAE_EXIT(CALCVELPOSFORNODES);
    },
    dash::tasks::in( this->deltatime()),
    dash::tasks::in( this->xdd(0)),
    dash::tasks::out(this->xd(0)),
    dash::tasks::out(this->x(0))
  );
  /*
  dash::tasks::ASYNC(
    [=](){
      const Real_t delt = deltatime();
      CalcPositionForNodes(delt, numNode());
    },
    dash::tasks::in(&this->deltatime()),
    dash::tasks::in(&this->xd(0)),
    dash::tasks::out(&this->x(0))
  );
  */

  // create communication tasks
  m_comm.Send_PosVel();
  m_comm.Sync_PosVel();
}

void Domain::LagrangeElements()
{
  // new relative vol -- temp
  dash::tasks::ASYNC(
    [&](){
      EXTRAE_ENTER(LAGRANGEELEMS);
      CalcLagrangeElements(m_vnew.data());
      EXTRAE_EXIT(LAGRANGEELEMS);
    },
    dash::tasks::in(this->x(0)),
    dash::tasks::out(this->dxx(0)),
    dash::tasks::out(this->vdov(0)),
    dash::tasks::out(this->arealg(0)),
    dash::tasks::out(m_vnew)
  );

  // Calculate Q.  (Monotonic q option requires communication)
  CalcQForElems(m_vnew);

  dash::tasks::ASYNC(
    [&](){
      EXTRAE_ENTER(MATERIALPROPERTIES);
      ApplyMaterialPropertiesForElems(m_vnew.data());

      UpdateVolumesForElems(m_vnew.data(), v_cut(), numElem());

      EXTRAE_EXIT(MATERIALPROPERTIES);
    },
    dash::tasks::in(this->x(0)),
    dash::tasks::out(this->v(0)),
    dash::tasks::out(this->p(0)),
    dash::tasks::out(this->e(0)),
    dash::tasks::out(this->q(0)),
    dash::tasks::out(this->ss(0)),
    dash::tasks::out(m_vnew)
  );
}


void Domain::CalcAccelerationForNodes(Index_t numNode)
{
//#pragma omp parallel for firstprivate(numNode)
      dash::tasks::TASKLOOP(0, numNode,
                            dash::tasks::num_chunks(dash::tasks::numthreads()*DASH_TASKLOOP_FACTOR),
        [&](Index_t from, Index_t to){
          for (Index_t i = from; i < to; ++i) {
            xdd(i) = fx(i) / nodalMass(i);
            ydd(i) = fy(i) / nodalMass(i);
            zdd(i) = fz(i) / nodalMass(i);
          }
        });
      dash::tasks::complete();
}

void Domain::ApplyAccelerationBoundaryConditionsForNodes()
{
  Domain &domain = (*this);
  Index_t size = domain.sizeX();
  Index_t numNodeBC = (size+1)*(size+1) ;

  auto grainsize = dash::tasks::num_chunks(dash::tasks::numthreads());

//#pragma omp parallel
  {
    if (!(domain.symmXempty() != 0) &&
        !(domain.symmYempty() != 0) &&
        !(domain.symmZempty() != 0))
    dash::tasks::TASKLOOP(Index_t{0}, numNodeBC, grainsize,
      [&](Index_t from, Index_t to){
        if (!domain.symmXempty() != 0)
          for(Index_t i=from; i<to; ++i)
            domain.xdd(domain.symmX(i)) = Real_t(0.0) ;

        if (!domain.symmYempty() != 0)
          for(Index_t i=from; i<to; ++i)
            domain.ydd(domain.symmY(i)) = Real_t(0.0) ;

        if (!domain.symmZempty() != 0)
          for(Index_t i=from; i<to; ++i)
            domain.zdd(domain.symmZ(i)) = Real_t(0.0) ;
      });
#if 0
    if (!domain.symmXempty() != 0) {
//#pragma omp for nowait firstprivate(numNodeBC)
      dash::tasks::TASKLOOP(Index_t{0}, numNodeBC, grainsize,
        [&](Index_t from, Index_t to){
          for(Index_t i=from; i<to; ++i)
            domain.xdd(domain.symmX(i)) = Real_t(0.0) ;
        });
    }

    if (!domain.symmYempty() != 0) {
//#pragma omp for nowait firstprivate(numNodeBC)
      dash::tasks::TASKLOOP(Index_t{0}, numNodeBC, grainsize,
        [&](Index_t from, Index_t to){
          for(Index_t i=from; i<to; ++i)
            domain.ydd(domain.symmY(i)) = Real_t(0.0) ;
        });
    }

    if (!domain.symmZempty() != 0) {
//#pragma omp for nowait firstprivate(numNodeBC)
      dash::tasks::TASKLOOP(Index_t{0}, numNodeBC, grainsize,
        [&](Index_t from, Index_t to){
          for(Index_t i=from; i<to; ++i)
            domain.zdd(domain.symmZ(i)) = Real_t(0.0) ;
        });
    }
#endif
  }
  dash::tasks::complete();
}

void Domain::CalcVelocityForNodes(const Real_t dt,
				  const Real_t u_cut,
				  Index_t numNode)
{
//#pragma omp parallel for firstprivate(numNode)
  dash::tasks::TASKLOOP(Index_t{0}, numNode,
                        dash::tasks::num_chunks(dash::tasks::numthreads()*DASH_TASKLOOP_FACTOR),
    [&, u_cut, dt](Index_t from, Index_t to){
      for ( Index_t i = from; i < to; ++i )
      {
        Real_t xdtmp, ydtmp, zdtmp ;

        xdtmp = xd(i) + xdd(i) * dt ;
        if( FABS(xdtmp) < u_cut ) xdtmp = Real_t(0.0);
        xd(i) = xdtmp ;

        ydtmp = yd(i) + ydd(i) * dt ;
        if( FABS(ydtmp) < u_cut ) ydtmp = Real_t(0.0);
        yd(i) = ydtmp ;

        zdtmp = zd(i) + zdd(i) * dt ;
        if( FABS(zdtmp) < u_cut ) zdtmp = Real_t(0.0);
        zd(i) = zdtmp ;
      }
    },
    [&](Index_t from, Index_t to, dash::tasks::DependencyVectorInserter deps){
      *deps = dash::tasks::out(xd(from));
      *deps = dash::tasks::in( xdd(from));
    });
  //dash::tasks::complete();
}


void Domain::CalcPositionForNodes(const Real_t dt,
				  Index_t numNode)
{
//#pragma omp parallel for firstprivate(numNode)
  dash::tasks::TASKLOOP(Index_t{0}, numNode,
                        dash::tasks::num_chunks(dash::tasks::numthreads()*DASH_TASKLOOP_FACTOR),
    [&, dt](Index_t from, Index_t to){
      for( Index_t i = from; i < to; ++i )
      {
        x(i) += xd(i) * dt;
        y(i) += yd(i) * dt;
        z(i) += zd(i) * dt;
      }
    },
    [&](Index_t from, Index_t to, dash::tasks::DependencyVectorInserter deps){
      *deps = dash::tasks::out(x(from));
      *deps = dash::tasks::in( xd(from));
    });
}


void Domain::CalcForceForNodes()
{

//#pragma omp parallel for firstprivate(numNode)
  dash::tasks::ASYNC(
    [&](){
      EXTRAE_ENTER(CALCFORCEFORNODES);
      Domain& domain = (*this);
      Index_t numNode = domain.numNode() ;
      dash::tasks::TASKLOOP(Index_t{0}, numNode,
                            dash::tasks::num_chunks(dash::tasks::numthreads()),
        [&](Index_t from, Index_t to){
          for (Index_t i=from; i<to; ++i) {
            domain.fx(i) = Real_t(0.0) ;
            domain.fy(i) = Real_t(0.0) ;
            domain.fz(i) = Real_t(0.0) ;
          }
        });
      dash::tasks::complete();

      // Calcforce calls partial, force, hourq
      CalcVolumeForceForElems(domain);
      EXTRAE_EXIT(CALCFORCEFORNODES);
    },
    dash::tasks::in(this->v(0)),
    dash::tasks::in(this->p(0)),
    dash::tasks::in(this->q(0)),
    dash::tasks::out(this->fx(0))
  );

  m_comm.Send_Force();
  m_comm.Sync_Force();
}

void Domain::AllocateStrains(Int_t numElem)
{
  extrae_task_event e(ALLOCATESTRAINS);
  m_dxx.resize(numElem);
  m_dyy.resize(numElem);
  m_dzz.resize(numElem);
}

void Domain::DeallocateStrains()
{
  m_dxx.clear();
  m_dyy.clear();
  m_dzz.clear();
}

void Domain::AllocateGradients(Int_t numElem, Int_t allElem)
{
  extrae_task_event e(ALLOCATEGRADIENTS);
  // Position gradients
  m_delx_xi.resize(numElem);
  m_delx_eta.resize(numElem) ;
  m_delx_zeta.resize(numElem) ;

  // Velocity gradients
  m_delv_xi.resize(allElem) ;
  m_delv_eta.resize(allElem);
  m_delv_zeta.resize(allElem) ;
}

void Domain::DeallocateGradients()
{
  m_delx_zeta.clear();
  m_delx_eta.clear();
  m_delx_xi.clear();

  m_delv_zeta.clear();
  m_delv_eta.clear();
  m_delv_xi.clear();
}


void Domain::SetupSymmetryPlanes(Int_t edgeNodes)
{
  //
  // for a process location on the border of the process grid
  // (plane==0 or col==0 or row==0 ) this routine computes the
  // node indices in z, y, x coordinate of the outermost plane
  //
  // NOTE: assumes a cube
  //

  // Boundary nodesets
  if (colLoc() == 0)
    m_symmX.resize(edgeNodes*edgeNodes);
  if (rowLoc() == 0)
    m_symmY.resize(edgeNodes*edgeNodes);
  if (planeLoc() == 0)
    m_symmZ.resize(edgeNodes*edgeNodes);

  Index_t nidx = 0 ;
  for (Index_t i=0; i<edgeNodes; ++i) {
    Index_t planeInc = i*edgeNodes*edgeNodes ;
    Index_t rowInc   = i*edgeNodes ;
    for (Index_t j=0; j<edgeNodes; ++j) {
      if (planeLoc() == 0) {
        m_symmZ[nidx] = rowInc   + j ;
      }
      if (rowLoc() == 0) {
        m_symmY[nidx] = planeInc + j ;
      }
      if (colLoc() == 0) {
        m_symmX[nidx] = planeInc + j*edgeNodes ;
      }
      ++nidx ;
    }
  }
}


void Domain::SetupThreadSupportStructures()
{
  Index_t numthreads = dash::tasks::numthreads();

  if (numthreads > 1) {
    // set up node-centered indexing of elements
    Index_t *nodeElemCount = new Index_t[numNode()] ;

    for (Index_t i=0; i<numNode(); ++i) {
      nodeElemCount[i] = 0 ;
    }

    for (Index_t i=0; i<numElem(); ++i) {
      Index_t *nl = nodelist(i) ;
      for (Index_t j=0; j < 8; ++j) {
        ++(nodeElemCount[nl[j]] );
      }
    }

    m_nodeElemStart = new Index_t[numNode()+1] ;

    m_nodeElemStart[0] = 0;
    for (Index_t i=1; i <= numNode(); ++i) {
      m_nodeElemStart[i] =
      m_nodeElemStart[i-1] + nodeElemCount[i-1] ;
    }

    m_nodeElemCornerList = new Index_t[m_nodeElemStart[numNode()]];

    for (Index_t i=0; i < numNode(); ++i) {
      nodeElemCount[i] = 0;
    }

    for (Index_t i=0; i < numElem(); ++i) {
      Index_t *nl = nodelist(i) ;
      for (Index_t j=0; j < 8; ++j) {
        Index_t m = nl[j];
        Index_t k = i*8 + j ;
        Index_t offset = m_nodeElemStart[m] + nodeElemCount[m] ;
        m_nodeElemCornerList[offset] = k;
        ++(nodeElemCount[m]) ;
      }
    }

    Index_t clSize = m_nodeElemStart[numNode()] ;
    for (Index_t i=0; i < clSize; ++i) {
      Index_t clv = m_nodeElemCornerList[i] ;
      if ((clv < 0) || (clv > numElem()*8)) {
	fprintf(stderr,
		"AllocateNodeElemIndexes(): nodeElemCornerList entry out of range!\n");
#if USE_MPI
	MPI_Abort(MPI_COMM_WORLD, -1);
#else
	exit(-1);
#endif
      }
    }

    delete [] nodeElemCount ;
  }
  else {
    // These arrays are not used if we're not threaded
    m_nodeElemStart = NULL;
    m_nodeElemCornerList = NULL;
  }
}

void Domain::InitializeFieldData()
{
  // Basic Field Initialization
  for (Index_t i=0; i<numElem(); ++i) {
    e(i) =  Real_t(0.0) ;
    p(i) =  Real_t(0.0) ;
    q(i) =  Real_t(0.0) ;
    ss(i) = Real_t(0.0) ;
  }

  // Note - v initializes to 1.0, not 0.0!
  for (Index_t i=0; i<numElem(); ++i) {
    v(i) = Real_t(1.0) ;
  }

  for (Index_t i=0; i<numNode(); ++i) {
    xd(i) = Real_t(0.0) ;
    yd(i) = Real_t(0.0) ;
    zd(i) = Real_t(0.0) ;
  }

  for (Index_t i=0; i<numNode(); ++i) {
    xdd(i) = Real_t(0.0) ;
    ydd(i) = Real_t(0.0) ;
    zdd(i) = Real_t(0.0) ;
  }

  for (Index_t i=0; i<numNode(); ++i) {
    nodalMass(i) = Real_t(0.0) ;
  }

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


void Domain::SetupElementConnectivities(Int_t edgeElems)
{
  m_lxim.resize(numElem());
  m_lxip.resize(numElem());
  m_letam.resize(numElem());
  m_letap.resize(numElem());
  m_lzetam.resize(numElem());
  m_lzetap.resize(numElem());

  lxim(0) = 0 ;
  for (Index_t i=1; i<numElem(); ++i) {
    lxim(i)   = i-1 ;
    lxip(i-1) = i ;
  }
  lxip(numElem()-1) = numElem()-1 ;

  for (Index_t i=0; i<edgeElems; ++i) {
    letam(i) = i ;
    letap(numElem()-edgeElems+i) = numElem()-edgeElems+i ;
  }
  for (Index_t i=edgeElems; i<numElem(); ++i) {
    letam(i) = i-edgeElems ;
    letap(i-edgeElems) = i ;
  }
  for (Index_t i=0; i<edgeElems*edgeElems; ++i) {
    lzetam(i) = i ;
    lzetap(numElem()-edgeElems*edgeElems+i) = numElem()-edgeElems*edgeElems+i ;
  }
  for (Index_t i=edgeElems*edgeElems; i<numElem(); ++i) {
    lzetam(i) = i - edgeElems*edgeElems ;
    lzetap(i-edgeElems*edgeElems) = i ;
  }
}


void Domain::SetupBoundaryConditions(Int_t edgeElems)
{
  Index_t ghostIdx[6] ;  // offsets to ghost locations

  Index_t rowMin   = (rowLoc()   == 0)       ? 0 : 1;
  Index_t rowMax   = (rowLoc()   == tp()-1)  ? 0 : 1;
  Index_t colMin   = (colLoc()   == 0)       ? 0 : 1;
  Index_t colMax   = (colLoc()   == tp()-1)  ? 0 : 1;
  Index_t planeMin = (planeLoc() == 0)       ? 0 : 1;
  Index_t planeMax = (planeLoc() == tp()-1)  ? 0 : 1;

  m_elemBC.resize(numElem());

  // set up boundary condition information
  for (Index_t i=0; i<numElem(); ++i) {
    elemBC(i) = Int_t(0) ;
  }

  for (Index_t i=0; i<6; ++i) {
    ghostIdx[i] = INT_MIN ;
  }

  Int_t pidx = numElem() ;
  if (planeMin != 0) {
    ghostIdx[0] = pidx ;
    pidx += sizeX()*sizeY() ;
  }

  if (planeMax != 0) {
    ghostIdx[1] = pidx ;
    pidx += sizeX()*sizeY() ;
  }

  if (rowMin != 0) {
    ghostIdx[2] = pidx ;
    pidx += sizeX()*sizeZ() ;
  }

  if (rowMax != 0) {
    ghostIdx[3] = pidx ;
    pidx += sizeX()*sizeZ() ;
  }

  if (colMin != 0) {
    ghostIdx[4] = pidx ;
    pidx += sizeY()*sizeZ() ;
  }

  if (colMax != 0) {
    ghostIdx[5] = pidx ;
  }

  // symmetry plane or free surface BCs
  for (Index_t i=0; i<edgeElems; ++i) {
    Index_t planeInc = i*edgeElems*edgeElems ;
    Index_t rowInc   = i*edgeElems ;
    for (Index_t j=0; j<edgeElems; ++j) {
      if (planeLoc() == 0) {
	elemBC(rowInc+j) |= ZETA_M_SYMM ;
      }
      else {
	elemBC(rowInc+j) |= ZETA_M_COMM ;
	lzetam(rowInc+j) = ghostIdx[0] + rowInc + j ;
      }

      if (planeLoc() == tp()-1) {
	elemBC(rowInc+j+numElem()-edgeElems*edgeElems) |=
	  ZETA_P_FREE;
      }
      else {
	elemBC(rowInc+j+numElem()-edgeElems*edgeElems) |=
	  ZETA_P_COMM ;
	lzetap(rowInc+j+numElem()-edgeElems*edgeElems) =
	  ghostIdx[1] + rowInc + j ;
      }

      if (rowLoc() == 0) {
	elemBC(planeInc+j) |= ETA_M_SYMM ;
      }
      else {
	elemBC(planeInc+j) |= ETA_M_COMM ;
	letam(planeInc+j) = ghostIdx[2] + rowInc + j ;
      }

      if (rowLoc() == tp()-1) {
	elemBC(planeInc+j+edgeElems*edgeElems-edgeElems) |=
	  ETA_P_FREE ;
      }
      else {
	elemBC(planeInc+j+edgeElems*edgeElems-edgeElems) |=
	  ETA_P_COMM ;
	letap(planeInc+j+edgeElems*edgeElems-edgeElems) =
	  ghostIdx[3] +  rowInc + j ;
      }

      if (colLoc() == 0) {
	elemBC(planeInc+j*edgeElems) |= XI_M_SYMM ;
      }
      else {
	elemBC(planeInc+j*edgeElems) |= XI_M_COMM ;
	lxim(planeInc+j*edgeElems) = ghostIdx[4] + rowInc + j ;
      }

      if (colLoc() == tp()-1) {
	elemBC(planeInc+j*edgeElems+edgeElems-1) |= XI_P_FREE ;
      }
      else {
	elemBC(planeInc+j*edgeElems+edgeElems-1) |= XI_P_COMM ;
	lxip(planeInc+j*edgeElems+edgeElems-1) =
	  ghostIdx[5] + rowInc + j ;
      }
    }
  }
}


void Domain::DepositInitialEnergy(Int_t nx)
{
  // deposit initial energy
  // An energy of 3.948746e+7 is correct for a problem with
  // 45 zones along a side - we need to scale it
  const Real_t ebase = Real_t(3.948746e+7);
  Real_t scale = (nx*tp())/Real_t(45.0);
  Real_t einit = ebase*scale*scale*scale;

  if (rowLoc() + colLoc() + planeLoc() == 0) {
    // Dump into the first zone (which we know is in the corner)
    // of the domain that sits at the origin
    e(0) = einit;
  }

  //set initial deltatime base on analytic CFL calculation
  deltatime() = (Real_t(.5)*cbrt(volo(0)))/sqrt(Real_t(2.0)*einit);
}


void Domain::CalcLagrangeElements(Real_t* vnew)
{
  Domain& domain = (*this);

  Index_t numElem = domain.numElem() ;
  if (numElem > 0) {
    const Real_t deltatime = domain.deltatime() ;

    //domain.AllocateStrains(numElem);

    CalcKinematicsForElems(domain, vnew, deltatime, numElem) ;

    // element loop to do some stuff not included in the elemlib function.
//#pragma omp parallel for firstprivate(numElem)

    dash::tasks::TASKLOOP(Index_t{0}, numElem,
                          dash::tasks::num_chunks(dash::tasks::numthreads()*DASH_TASKLOOP_FACTOR),
    [&](Index_t from, Index_t to){
      for ( Index_t k=from; k<to; ++k )
      {
        // calc strain rate and apply as constraint (only done in FB element)
        Real_t vdov = domain.dxx(k) + domain.dyy(k) + domain.dzz(k) ;
        Real_t vdovthird = vdov/Real_t(3.0) ;

        // make the rate of deformation tensor deviatoric
        domain.vdov(k) = vdov ;
        domain.dxx(k) -= vdovthird ;
        domain.dyy(k) -= vdovthird ;
        domain.dzz(k) -= vdovthird ;

        // See if any volumes are negative, and take appropriate action.
        if (vnew[k] <= Real_t(0.0))
        {
          std::cerr << "VOLUME ERROR: vnew[k]=" << vnew[k] << std::endl;
#if USE_MPI
          MPI_Abort(MPI_COMM_WORLD, VolumeError) ;
#elif USE_DASH
          dart_abort(VolumeError);
#else
          exit(VolumeError);
#endif
        }
      }
      });
    dash::tasks::complete();
    //domain.DeallocateStrains();
  }
}

void Domain::CalcQForElems(std::vector<Real_t>& vnew)
{
  //
  // MONOTONIC Q option
  //

  Domain& domain = (*this);

  if (domain.numElem() != 0) {
    dash::tasks::ASYNC(
      [=, &domain, &vnew](){
        EXTRAE_ENTER(MONOQGRADIENTS);

#if 0
        Index_t numElem = domain.numElem() ;
        Int_t allElem = numElem +           /* local elem */
          2*domain.sizeX()*domain.sizeY() + /* plane ghosts */
          2*domain.sizeX()*domain.sizeZ() + /* row ghosts */
          2*domain.sizeY()*domain.sizeZ() ; /* col ghosts */
        domain.AllocateGradients(numElem, allElem);
#endif

        /* Calculate velocity gradients */
        CalcMonotonicQGradientsForElems(domain, vnew.data());
        EXTRAE_EXIT(MONOQGRADIENTS);
      },
      dash::tasks::in(vnew),
      dash::tasks::out(this->delv_xi(0))
    );

    m_comm.Send_MonoQ();
    m_comm.Sync_MonoQ();

    dash::tasks::ASYNC(
      [=, &domain, &vnew](){
        EXTRAE_ENTER(MONOQELEMS);
        CalcMonotonicQForElems(domain, vnew.data()) ;

        // Free up memory
        //domain.DeallocateGradients();

        /* Don't allow excessive artificial viscosity */
        Index_t idx = -1;
        Index_t numElem = domain.numElem() ;
        for (Index_t i=0; i<numElem; ++i) {
          if ( domain.q(i) > domain.qstop() ) {
            idx = i ;
            break ;
          }
        }

        if(idx >= 0) {
          std::cerr << "VOLUME ERROR: domain.q(" << idx << ")=" << domain.q(idx)
                    << std::endl;
#if USE_MPI
          MPI_Abort(MPI_COMM_WORLD, QStopError) ;
#elif USE_DASH
          dart_abort(VolumeError);
#else
          exit(QStopError);
#endif
        }
        EXTRAE_EXIT(MONOQELEMS);
      },
      dash::tasks::in(this->delv_xi(0)),
      dash::tasks::out(vnew),
      dash::tasks::out(domain.qq(0)),
      dash::tasks::out(domain.ql(0))
    );
  }
}

void Domain::ApplyMaterialPropertiesForElems(Real_t vnew[])
{
  Domain& domain = (*this);
  Index_t numElem = domain.numElem() ;

  if (numElem != 0) {
    /* Expose all of the variables needed for material evaluation */
    Real_t eosvmin = domain.eosvmin() ;
    Real_t eosvmax = domain.eosvmax() ;

//#pragma omp for firstprivate(numElem)
    dash::tasks::TASKLOOP(Index_t{0}, numElem,
                          dash::tasks::num_chunks(dash::tasks::numthreads()),
      [&, eosvmin](Index_t from, Index_t to){
        if (eosvmin != Real_t(0.))
          for(Index_t i=from; i<to; ++i) {
            if (vnew[i] < eosvmin)
              vnew[i] = eosvmin ;
          }

        if (eosvmax != Real_t(0.))
          for(Index_t i=from; i<to; ++i) {
            if (vnew[i] > eosvmax)
              vnew[i] = eosvmax ;
          }
      });
    //dash::tasks::complete();

//#pragma omp parallel
    {
#if 0
      // Bound the updated relative volumes with eosvmin/max
      if (eosvmin != Real_t(0.)) {
//#pragma omp for firstprivate(numElem)
        dash::tasks::TASKLOOP(Index_t{0}, numElem,
          [&, eosvmin](Index_t from, Index_t to){
            for(Index_t i=from; i<to; ++i) {
              if (vnew[i] < eosvmin)
                vnew[i] = eosvmin ;
            }
          },
          [vnew]
          (Index_t from, Index_t to, dash::tasks::DependencyVectorInserter deps){
            *deps = dash::tasks::out(&vnew[from]);
          });
        //dash::tasks::complete();
      }

      if (eosvmax != Real_t(0.)) {
//#pragma omp for nowait firstprivate(numElem)
        dash::tasks::TASKLOOP(Index_t{0}, numElem,
                              numElem/(dash::tasks::numthreads()*DASH_TASKLOOP_FACTOR),
          [&, eosvmax](Index_t from, Index_t to){
            for(Index_t i=from; i<to; ++i) {
              if (vnew[i] > eosvmax)
                vnew[i] = eosvmax ;
            }
          },
          [vnew]
          (Index_t from, Index_t to, dash::tasks::DependencyVectorInserter deps){
            *deps = dash::tasks::out(&vnew[from]);
          });
        //dash::tasks::complete();
      }
#endif

      // This check may not make perfect sense in LULESH, but
      // it's representative of something in the full code -
      // just leave it in, please
//#pragma omp for nowait firstprivate(numElem)
      dash::tasks::TASKLOOP(Index_t{0}, numElem,
                            dash::tasks::num_chunks(dash::tasks::numthreads()),
        [&, eosvmin, eosvmax](Index_t from, Index_t to){
          for (Index_t i=from; i<to; ++i) {
            Real_t vc = domain.v(i) ;
            if (eosvmin != Real_t(0.)) {
              if (vc < eosvmin)
                vc = eosvmin ;
            }
            if (eosvmax != Real_t(0.)) {
              if (vc > eosvmax)
                vc = eosvmax ;
            }
            if (vc <= 0.) {
              std::cerr << "VOLUME ERROR: vc=" << vc << std::endl;
#if USE_MPI
              MPI_Abort(MPI_COMM_WORLD, VolumeError) ;
#elif USE_DASH
              dart_abort(VolumeError);
#else
              exit(VolumeError);
#endif
            }
          }
        });
      dash::tasks::complete();
    }

    dash::tasks::ASYNC(
      [=, &domain](){
        for (Int_t r=0 ; r<domain.numReg() ; r++) {
          Index_t numElemReg = domain.regElemSize(r);
          Index_t *regElemList = domain.regElemlist(r);
          Int_t rep;
          //Determine load imbalance for this region
          //round down the number with lowest cost
          if(r < domain.numReg()/2)
            rep = 1;
          //you don't get an expensive region unless you at least have 5 regions
          else if(r < (domain.numReg() - (domain.numReg()+15)/20))
            rep = 1 + domain.cost();
          //very expensive regions
          else
            rep = 10 * (1+ domain.cost());
          EvalEOSForElems(domain, vnew, numElemReg, regElemList, rep);
        }
      }
    );
    dash::tasks::complete();
  }
}

void Domain::UpdateVolumesForElems(Real_t *vnew,
			   Real_t v_cut, Index_t length)
{
  Domain& domain = (*this);

  if (length != 0) {
//#pragma omp parallel for firstprivate(length, v_cut)
    dash::tasks::TASKLOOP(Index_t{0}, length,
                          dash::tasks::num_chunks(dash::tasks::numthreads()),
      [&, v_cut, vnew](Index_t from, Index_t to){
        for(Index_t i=from; i<to; ++i) {
          Real_t tmpV = vnew[i] ;

          if ( FABS(tmpV - Real_t(1.0)) < v_cut )
            tmpV = Real_t(1.0) ;

          domain.v(i) = tmpV ;
        }
      });
    dash::tasks::complete();
  }
}



void Domain::VerifyAndWriteFinalOutput(Real_t elapsed,
				       Int_t  nx,
				       Int_t  numRanks)
{
  // GrindTime1 only takes a single domain into account, and is thus a
  // good way to measure processor speed indepdendent of MPI
  // parallelism.
  //
  // GrindTime2 takes into account speedups from MPI parallelism
  Real_t grindTime1 = ((elapsed*1.0e6)/cycle())/(nx*nx*nx);
  Real_t grindTime2 = ((elapsed*1.0e6)/cycle())/(nx*nx*nx*numRanks);

  Index_t ElemId = 0;

#if 0
  cout << "Run completed:" << endl;
  cout << "   Problem size        = " << nx << endl;
  cout << "   MPI tasks           = " << numRanks << endl;
  cout << "   Iteration count     = " << cycle() << endl;
  cout << "   Final Origin Energy = " << e(ElemId) << endl;
#endif

  printf("Run completed:  \n");
  printf("   Problem size        =  %i \n",    nx);
  printf("   MPI tasks           =  %i \n",    numRanks);
  printf("   Iteration count     =  %i \n",    cycle());
  printf("   Final Origin Energy = %12.6e \n", e(ElemId));

  Real_t   MaxAbsDiff = Real_t(0.0);
  Real_t TotalAbsDiff = Real_t(0.0);
  Real_t   MaxRelDiff = Real_t(0.0);

  for (Index_t j=0; j<nx; ++j) {
    for (Index_t k=j+1; k<nx; ++k) {
      Real_t AbsDiff = FABS(e(j*nx+k)-e(k*nx+j));
      TotalAbsDiff  += AbsDiff;

      if (MaxAbsDiff <AbsDiff) MaxAbsDiff = AbsDiff;

      Real_t RelDiff = AbsDiff / e(k*nx+j);

      if (MaxRelDiff <RelDiff)  MaxRelDiff = RelDiff;
    }
  }

#if 0
  cout << std::scientific;
  // Quick symmetry check
  cout << "   Testing Plane 0 of Energy Array on rank 0: "<< endl;
  cout << "        MaxAbsDiff   = " << MaxAbsDiff   << endl;
  cout << "        TotalAbsDiff = " << TotalAbsDiff << endl;
  cout << "        MaxRelDiff   = " << MaxRelDiff   << endl;

  // Timing information
  cout << endl;
  cout << "Elapsed time         = " << elapsed << endl;
  cout << "Grind time (us/z/c)  = " << grindTime1 << " (per dom) " << endl;
  cout << "Grind time (us/z/c)  = " << grindTime2 << " (overall) " << endl;
  cout << "FOM                  = " << 1000.0/grindTime2 << endl;
#endif
  printf("   Testing Plane 0 of Energy Array on rank 0:\n");
  printf("        MaxAbsDiff   = %12.6e\n",   MaxAbsDiff   );
  printf("        TotalAbsDiff = %12.6e\n",   TotalAbsDiff );
  printf("        MaxRelDiff   = %12.6e\n\n", MaxRelDiff   );

   // Timing information
  printf("\nElapsed time         = %10.2f (s)\n", elapsed);
  printf("Grind time (us/z/c)  = %10.8g (per dom)  (%10.8g overall)\n", grindTime1, grindTime2);
  printf("FOM                  = %10.8g (z/s)\n\n", 1000.0/grindTime2); // zones per second
}
