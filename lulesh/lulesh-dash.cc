
#include <libdash.h>
#include <iostream>
#include "lulesh-dash.h"
#include "lulesh-calc.h"
#ifdef USE_MPI
#include "lulesh-mpi.h"
#endif

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


RegionIndexSet::RegionIndexSet(Int_t nr, Int_t balance, Index_t numElem)
{
#if USE_MPI
  Index_t myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;
  srand(myRank);
#else
  srand(0);
  Index_t myRank = 0;
#endif

  this->numReg() = nr;
  m_regElemSize = new Index_t[numReg()];
  m_regElemlist = new Index_t*[numReg()];
  m_regNumList  = new Index_t[numElem] ;  // material indexset

  Index_t nextIndex = 0;

  //if we only have one region just fill it
  // Fill out the regNumList with material numbers, which are always
  // the region index plus one
  if(numReg() == 1) {
    while (nextIndex < numElem) {
      this->regNumList(nextIndex) = 1;
      nextIndex++;
    }
    regElemSize(0) = 0;
  }
  //If we have more than one region distribute the elements.
  else {
    Int_t regionNum;
    Int_t regionVar;
    Int_t lastReg = -1;
    Int_t binSize;
    Index_t elements;
    Index_t runto = 0;
    Int_t costDenominator = 0;
    Int_t* regBinEnd = new Int_t[numReg()];
    //Determine the relative weights of all the regions.  This is
    //based off the -b flag.  Balance is the value passed into b.
    for (Index_t i=0 ; i<numReg() ; ++i) {
      regElemSize(i) = 0;
      costDenominator += pow((i+1), balance);  //Total sum of all
					       //regions weights
      regBinEnd[i] = costDenominator;  //Chance of hitting a given
				       //region is (regBinEnd[i] -
				       //regBinEdn[i-1])/costDenominator
    }

    //Until all elements are assigned
    while (nextIndex < numElem) {
      //pick the region
      regionVar = rand() % costDenominator;
      Index_t i = 0;
      while(regionVar >= regBinEnd[i])
	i++;

      //rotate the regions based on MPI rank.  Rotation is Rank %
      //NumRegions this makes each domain have a different region with
      //the highest representation
      regionNum = ((i + myRank) % numReg()) + 1;
      // make sure we don't pick the same region twice in a row
      while(regionNum == lastReg) {
	regionVar = rand() % costDenominator;
	i = 0;
	while(regionVar >= regBinEnd[i])
	  i++;
	regionNum = ((i + myRank) % numReg()) + 1;
      }
      //Pick the bin size of the region and determine the number of
      //elements.
      binSize = rand() % 1000;
      if(binSize < 773) {
	elements = rand() % 15 + 1;
      }
      else if(binSize < 937) {
	elements = rand() % 16 + 16;
      }
      else if(binSize < 970) {
	elements = rand() % 32 + 32;
      }
      else if(binSize < 974) {
	elements = rand() % 64 + 64;
      }
      else if(binSize < 978) {
	elements = rand() % 128 + 128;
      }
      else if(binSize < 981) {
	elements = rand() % 256 + 256;
      }
      else
	elements = rand() % 1537 + 512;
      runto = elements + nextIndex;

      //Store the elements.  If we hit the end before we run out of elements then just stop.
      while (nextIndex < runto && nextIndex < numElem) {
	this->regNumList(nextIndex) = regionNum;
	nextIndex++;
      }
      lastReg = regionNum;
    }
  }

  // Convert regNumList to region index sets
  // First, count size of each region
  for (Index_t i=0 ; i<numElem ; ++i) {
    int r = this->regNumList(i)-1; // region index == regnum-1
    regElemSize(r)++;
  }

  // Second, allocate each region index set
  for (Index_t i=0 ; i<numReg() ; ++i) {
    m_regElemlist[i] = new Index_t[regElemSize(i)];
    regElemSize(i) = 0;
  }

  // Third, fill index sets
  for (Index_t i=0 ; i<numElem ; ++i) {
    Index_t r = regNumList(i)-1;       // region index == regnum-1
    Index_t regndx = regElemSize(r)++; // Note increment
    regElemlist(r,regndx) = i;
  }
}

RegionIndexSet::~RegionIndexSet()
{
#if 0
  for (Index_t i=0 ; i<numReg() ; ++i) {
    delete[](m_regElemlist[i]);
  }
  delete[](m_regElemSize);
  delete[](m_regElemlist);
#endif
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

  m_region(opts.numReg(), opts.balance(), numElem() ),

  m_comm(*this)
{
  //
  // node - centered fields
  //
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

  //
  // element-centered fields
  //
  m_e.allocate(m_ElemPat);
  m_p.allocate(m_ElemPat);
  m_q.allocate(m_ElemPat);
  m_ql.allocate(m_ElemPat);
  m_qq.allocate(m_ElemPat);

  m_v.allocate(m_ElemPat);
  m_volo.allocate(m_ElemPat);
  m_delv.allocate(m_ElemPat);
  m_vdov.allocate(m_ElemPat);

#if 0
  m_nodelist.allocate(m_ElemPat);
  m_volo.allocate(m_ElemPat);
  m_elemMass.allocate(m_ElemPat);

  m_dxx.resize(numElem());
  m_dyy.resize(numElem());
  m_dzz.resize(numElem());

  m_delv_xi.resize(numElem());
  m_delv_eta.resize(numElem());
  m_delv_zeta.resize(numElem());

  m_delx_xi.resize(numElem());
  m_delx_eta.resize(numElem());
  m_delx_zeta.resize(numElem());

  m_v.resize(numElem());
  m_vnew.resize(numElem());
  m_delv.resize(numElem());
  m_vdov.resize(numElem());

  m_arealg.resize(numElem());
  m_ss.resize(numElem());
#endif

  Index_t edgeNodes = m_nodelist.local.extent(0);

  BuildMesh();
  // XXX SetupThreadSupportStructures();

  // XXX CreateRegionIndexSets();

  // Setup symmetry nodesets
  SetupSymmetryPlanes(edgeNodes);

  // Setup element connectivities
  SetupElementConnectivities(edgeElems);

  // Setup symmetry planes and free surface boundary arrays
  SetupBoundaryConditions(edgeElems);

  InitializeFieldData();

  DepositInitialEnergy();

  /*
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
  */
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
  os << "[ " << myid << "] colLoc()          : " << colLoc()       << std::endl;
  os << "[ " << myid << "] rowLoc()          : " << rowLoc()       << std::endl;
  os << "[ " << myid << "] planeLoc()        : " << planeLoc()     << std::endl;
  os << "[ " << myid << "] numElem()         : " << numElem()      << std::endl;
  os << "[ " << myid << "] numElem(0)        : " << numElem(0)     << std::endl;
  os << "[ " << myid << "] numElem(1)        : " << numElem(1)     << std::endl;
  os << "[ " << myid << "] numElem(2)        : " << numElem(2)     << std::endl;
  os << "[ " << myid << "] numNode()         : " << numNode()      << std::endl;
  os << "[ " << myid << "] numNode(0)        : " << numNode(0)     << std::endl;
  os << "[ " << myid << "] numNode(1)        : " << numNode(1)     << std::endl;
  os << "[ " << myid << "] numNode(2)        : " << numNode(2)     << std::endl;
  os << "[ " << myid << "] sizeX()           : " << sizeX()        << std::endl;
  os << "[ " << myid << "] sizeY()           : " << sizeY()        << std::endl;
  os << "[ " << myid << "] sizeZ()           : " << sizeZ()        << std::endl;
  os << "[ " << myid << "] maxEdgeSize()     : " << maxEdgeSize()  << std::endl;
  os << "[ " << myid << "] maxPlaneSize()    : " << maxPlaneSize() << std::endl;

  os << "[ " << myid << "] time()            : " << time()            << std::endl;
  os << "[ " << myid << "] deltatime()       : " << deltatime()       << std::endl;
  os << "[ " << myid << "] deltatimemultlb() : " << deltatimemultlb() << std::endl;
  os << "[ " << myid << "] deltatimemultub() : " << deltatimemultub() << std::endl;
}


void Domain::InitialBoundaryExchange()
{
  m_comm.ExchangeNodalMass();
}


void Domain::TimeIncrement()
{
  Domain& domain = (*this);
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

#if USE_MPI
    MPI_Allreduce(&gnewdt, &newdt, 1,
		  ((sizeof(Real_t) == 4) ? MPI_FLOAT : MPI_DOUBLE),
		  MPI_MIN, MPI_COMM_WORLD) ;
#else
    newdt = gnewdt;
#endif

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
}


void Domain::LagrangeLeapFrog()
{
  LagrangeNodal();
  LagrangeElements();

#if SEDOV_SYNC_POS_VEL_LATE
  // initiate communication
  m_comm.Recv_PosVel();
  m_comm.Send_PosVel();
#endif

  CalcTimeConstraintsForElems(*this);

#if SEDOV_SYNC_POS_VEL_LATE
  // wait for completion
  m_comm.Sync_PosVel();
#endif
}

void Domain::LagrangeNodal()
{
  const Real_t delt = deltatime();
  Real_t ucut = u_cut();

  CalcForceForNodes();

#ifdef SEDOV_SYNC_POS_VEL_EARLY
  // post receives
  m_comm.Recv_PosVel();
#endif

  CalcAccelerationForNodes(numNode());
  ApplyAccelerationBoundaryConditionsForNodes();
  CalcVelocityForNodes(delt, ucut, numNode()) ;
  CalcPositionForNodes(delt, numNode());

#ifdef SEDOV_SYNC_POS_VEL_EARLY
  // finish communication
  m_comm.Send_PosVel();
  m_comm.Sync_PosVel();
#endif
}

void Domain::LagrangeElements()
{
  // new relative vol -- temp
  Real_t *vnew = Allocate<Real_t>(numElem());
  Domain &domain = (*this);

  /*
  CalcLagrangeElements(domain, vnew);

  // Calculate Q.  (Monotonic q option requires communication)
  CalcQForElems(domain, vnew);

  ApplyMaterialPropertiesForElems(domain, vnew);

  UpdateVolumesForElems(domain, vnew,
                        domain.v_cut(), numElem);
  */

  Release(&vnew);
}


void Domain::CalcAccelerationForNodes(Index_t numNode)
{
#pragma omp parallel for firstprivate(numNode)
  for (Index_t i = 0; i < numNode; ++i) {
    xdd(i) = fx(i) / nodalMass(i);
    ydd(i) = fy(i) / nodalMass(i);
    zdd(i) = fz(i) / nodalMass(i);
  }
}

void Domain::ApplyAccelerationBoundaryConditionsForNodes()
{
  Domain& domain    = (*this);
  Index_t size      = domain.sizeX();
  Index_t numNodeBC = (size+1)*(size+1) ;

#pragma omp parallel
  {
    if (!domain.symmXempty() != 0) {
#pragma omp for nowait firstprivate(numNodeBC)
      for(Index_t i=0 ; i<numNodeBC ; ++i)
	domain.xdd(domain.symmX(i)) = Real_t(0.0) ;
    }

    if (!domain.symmYempty() != 0) {
#pragma omp for nowait firstprivate(numNodeBC)
      for(Index_t i=0 ; i<numNodeBC ; ++i)
	domain.ydd(domain.symmY(i)) = Real_t(0.0) ;
    }

    if (!domain.symmZempty() != 0) {
#pragma omp for nowait firstprivate(numNodeBC)
      for(Index_t i=0 ; i<numNodeBC ; ++i)
	domain.zdd(domain.symmZ(i)) = Real_t(0.0) ;
    }
  }

}

void Domain::CalcVelocityForNodes(const Real_t dt,
				  const Real_t u_cut,
				  Index_t numNode)
{
#pragma omp parallel for firstprivate(numNode)
  for ( Index_t i = 0 ; i < numNode ; ++i )
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
}


void Domain::CalcPositionForNodes(const Real_t dt,
				  Index_t numNode)
{
#pragma omp parallel for firstprivate(numNode)
  for( Index_t i = 0 ; i < numNode ; ++i )
    {
      x(i) += xd(i) * dt;
      y(i) += yd(i) * dt;
      z(i) += zd(i) * dt;
    }
}


void Domain::CalcForceForNodes()
{
  Domain& domain = (*this);
  Index_t numNode = domain.numNode() ;

  m_comm.Recv_Force();

#pragma omp parallel for firstprivate(numNode)
  for (Index_t i=0; i<numNode; ++i) {
    domain.fx(i) = Real_t(0.0) ;
    domain.fy(i) = Real_t(0.0) ;
    domain.fz(i) = Real_t(0.0) ;
  }

  // Calcforce calls partial, force, hourq
  CalcVolumeForceForElems(domain) ;

  m_comm.Send_Force();
  m_comm.Sync_Force();
}


void Domain::AllocateStrains(Int_t numElem)
{
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
  // Position gradients
  m_delx_xi.resize(numElem) ;
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
