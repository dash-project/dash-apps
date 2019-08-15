#ifndef LULESH_DASH_H_INCLUDED
#define LULESH_DASH_H_INCLUDED

#include <libdash.h>
#include "lulesh.h"
#include "lulesh-opts.h"
#include "lulesh-dash-params.h"
#include "lulesh-dash-regions.h"
#ifdef USE_MPI
#include "lulesh-comm-mpi.h"
#endif
#ifdef USE_DASH
#include "lulesh-comm-dash.h"
#endif

/*
 * The DASH version of LULESH's 'Domain' data structure.
 *
 * In the MPI version this holds local data in std::vector containers,
 * here we use DASH Matrix. In the MPI version the domain is purely
 * local. Using DASH we have a global view of data from all processes,
 * but of course can also emulate the local view using the .local
 * accessor object.
 */
class Domain
{
private:
  using ElemPatternT = dash::Pattern<3>;
  using NodePatternT = dash::Pattern<3>;

  template<typename ElemType>
  using NodeMatrixT = std::vector<ElemType>;
  //using NodeMatrixT = dash::Matrix<ElemType, 3,
  //          NodePatternT::index_type,
  //          NodePatternT>;

  template<typename ElemType>
  using ElemMatrixT = std::vector<ElemType>;
  //using ElemMatrixT = dash::Matrix<ElemType, 3,
  //          ElemPatternT::index_type,
  //          ElemPatternT>;

  // all simulation constanst and parameters packed into a separate
  // struct for clarity
  Parameters m_param;

  // number of local elements and nodes in each dimension
  std::array<Index_t, 3> m_nElem;
  std::array<Index_t, 3> m_nNode;

  //
  // LULESH uses the terms 'col', 'row', and 'plane' to refer to the
  // 3D position of a process in the process grid. To get the same
  // ordering in DASH, a TeamSpec and the following equvalencies can
  // be used:
  //
  // m_ts.z() <=> col
  // m_ts.y() <=> row
  // m_ts.x() <=> plane
  dash::TeamSpec<3> m_ts;

  // pattern for element-centered data
  ElemPatternT m_ElemPat;

  // pattern for node-centered data
  NodePatternT m_NodePat;

  //
  // node-centered fields
  //
  NodeMatrixT<Real_t> m_x  , m_y  , m_z  ;   // coordinates
  NodeMatrixT<Real_t> m_xd , m_yd , m_zd ;   // velocities
  NodeMatrixT<Real_t> m_xdd, m_ydd, m_zdd;   // accelerations
  NodeMatrixT<Real_t> m_fx , m_fy , m_fz ;   // forces
  NodeMatrixT<Real_t> m_nodalMass;           // mass

  //
  // element-centered fields
  //
  ElemMatrixT<Real_t> m_e;         // energy
  ElemMatrixT<Real_t> m_p;         // pressure
  ElemMatrixT<Real_t> m_q;         // artificial viscosity
  ElemMatrixT<Real_t> m_ql;        // linear term for q
  ElemMatrixT<Real_t> m_qq;        // quadratic term for q

  ElemMatrixT<Real_t> m_v;         // relative volume
  ElemMatrixT<Real_t> m_volo;      // reference volume
  std::vector<Real_t> m_vnew ;     /* new relative volume -- temporary */

  ElemMatrixT<Real_t> m_delv;      // vnew - m_v (vnew is temp. allocated)
  ElemMatrixT<Real_t> m_vdov;      // volume derivative over volume
  ElemMatrixT<Real_t> m_elemMass;  // mass

  // elemToNode connectivity
  using NodeVecT = std::array<Index_t, 8>;
  ElemMatrixT<NodeVecT>  m_nodelist;

  // element connectivity across each face
  std::vector<Index_t>  m_lxim;
  std::vector<Index_t>  m_lxip;
  std::vector<Index_t>  m_letam;
  std::vector<Index_t>  m_letap;
  std::vector<Index_t>  m_lzetam;
  std::vector<Index_t>  m_lzetap;

  // symmetry/free-surface flags for each elem face
  std::vector<Int_t>    m_elemBC;

  Real_t               *m_dxx ;  /* principal strains -- temporary */
  Real_t               *m_dyy ;
  Real_t               *m_dzz ;

  Real_t               *m_delv_xi ;    /* velocity gradient -- temporary */
  Real_t               *m_delv_eta ;
  Real_t               *m_delv_zeta ;

  Real_t               *m_delx_xi ;    /* coordinate gradient -- temporary */
  Real_t               *m_delx_eta ;
  Real_t               *m_delx_zeta ;


  //
  // XXX - not really necessary for DASH
  //
  std::vector<Index_t> m_symmX;  // symmetry plane nodesets
  std::vector<Index_t> m_symmY;
  std::vector<Index_t> m_symmZ;

  RegionIndexSet m_region;

  std::vector<Real_t> m_arealg ;  // characteristic length of an element
  std::vector<Real_t> m_ss ;      // "sound speed"

  // OMP hack
  Index_t *m_nodeElemStart ;
  Index_t *m_nodeElemCornerList ;

#ifdef USE_MPI
  Comm m_comm;
#endif
#ifdef USE_DASH
  DASHComm m_comm;
#endif

  Int_t m_chunksize;

private:
  // helper routines used in constructor
  void BuildMesh();

  void SetupThreadSupportStructures();
  void SetupSymmetryPlanes(Int_t edgeNodes);
  void SetupElementConnectivities(Int_t edgeElems);
  void SetupBoundaryConditions(Int_t edgeElems);

  void CreateRegionIndexSets(Int_t nreg, Int_t balance);
  void InitializeFieldData();
  void DepositInitialEnergy(Int_t nx);

  void AllocateStrains(Int_t numElem);
  void DeallocateStrains();
  void AllocateVnew(Int_t numElem);
  void DeallocateVnew();
  void AllocateGradients(Int_t numElem, Int_t allElem);
  void DeallocateGradients();

public:
  void PrintDomain(){
    std::cout << "DOMAIN: \n  m_x:\n";
    for (auto ptr = m_x.begin(); ptr != m_x.end(); ptr++) {
      std::cout << *ptr << " ";
    }
    std::cout << std::endl;
    std::cout << " m_y:";
    for (auto ptr = m_y.begin(); ptr != m_y.end(); ptr++) {
      std::cout << *ptr << " ";
    }
    std::cout << std::endl;
    std::cout << " m_z:";
    for (auto ptr = m_z.begin(); ptr != m_z.end(); ptr++) {
      std::cout << *ptr << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }

public:
  Domain(const CmdLineOpts& opts);
  ~Domain();

  double wtime() { return m_comm.wtime(); }
  template<typename T> T allreduce_min(T val);
  template<typename T> T reduce_max(T val);

  // nodal coordinates
  Real_t& x(Index_t idx)     { return m_x[idx]; }
  Real_t& y(Index_t idx)     { return m_y[idx]; }
  Real_t& z(Index_t idx)     { return m_z[idx]; }

  // nodal velocities
  Real_t& xd(Index_t idx)    { return m_xd[idx]; }
  Real_t& yd(Index_t idx)    { return m_yd[idx]; }
  Real_t& zd(Index_t idx)    { return m_zd[idx]; }

  // nodal accelerations
  Real_t& xdd(Index_t idx)   { return m_xdd[idx]; }
  Real_t& ydd(Index_t idx)   { return m_ydd[idx]; }
  Real_t& zdd(Index_t idx)   { return m_zdd[idx]; }

  // nodal forces
  Real_t& fx(Index_t idx)    { return m_fx[idx]; }
  Real_t& fy(Index_t idx)    { return m_fy[idx]; }
  Real_t& fz(Index_t idx)    { return m_fz[idx]; }

  // nodal mass
  Real_t& nodalMass(Index_t idx)  { return m_nodalMass[idx]; }

  // element mass
  Real_t& elemMass(Index_t idx)   { return m_elemMass[idx]; }

  // energy, pressure, viscosity
  Real_t& e(Index_t idx)     { return m_e[idx]; }
  Real_t& p(Index_t idx)     { return m_p[idx]; }
  Real_t& q(Index_t idx)     { return m_q[idx]; }
  Real_t& ql(Index_t idx)    { return m_ql[idx]; }
  Real_t& qq(Index_t idx)    { return m_qq[idx]; }

  // volume, ref. volume, rel. volume, vol. derivative
  Real_t& v(Index_t idx)     { return m_v[idx]; }
  Real_t& volo(Index_t idx)  { return m_volo[idx]; }
  Real_t& delv(Index_t idx)  { return m_delv[idx]; }
  Real_t& vdov(Index_t idx)  { return m_vdov[idx]; }

  // principal strains -- temporary
  Real_t& dxx(Index_t idx)  { return m_dxx[idx] ; }
  Real_t& dyy(Index_t idx)  { return m_dyy[idx] ; }
  Real_t& dzz(Index_t idx)  { return m_dzz[idx] ; }

  // New relative volume - temporary
  Real_t& vnew(Index_t idx)  { return m_vnew[idx] ; }


  // velocity gradient -- temporary
  Real_t& delv_xi(Index_t idx)    { return m_delv_xi[idx] ; }
  Real_t& delv_eta(Index_t idx)   { return m_delv_eta[idx] ; }
  Real_t& delv_zeta(Index_t idx)  { return m_delv_zeta[idx] ; }

  // position gradient -- temporary
  Real_t& delx_xi(Index_t idx)    { return m_delx_xi[idx] ; }
  Real_t& delx_eta(Index_t idx)   { return m_delx_eta[idx] ; }
  Real_t& delx_zeta(Index_t idx)  { return m_delx_zeta[idx] ; }

  // elemToNode connectivity
  Index_t* nodelist(Index_t idx)  { return &m_nodelist[idx][0]; }

  // elem connectivities through face
  Index_t& lxim(  Index_t idx) { return m_lxim[idx]; }
  Index_t& lxip(  Index_t idx) { return m_lxip[idx]; }
  Index_t& letam( Index_t idx) { return m_letam[idx]; }
  Index_t& letap( Index_t idx) { return m_letap[idx]; }
  Index_t& lzetam(Index_t idx) { return m_lzetam[idx]; }
  Index_t& lzetap(Index_t idx) { return m_lzetap[idx]; }

  // elem face symm/free-surface flag
  Int_t&  elemBC(Index_t idx)  { return m_elemBC[idx]; }

  // nodes on symmetry planes
  Index_t symmX(Index_t idx) { return m_symmX[idx] ; }
  Index_t symmY(Index_t idx) { return m_symmY[idx] ; }
  Index_t symmZ(Index_t idx) { return m_symmZ[idx] ; }
  bool symmXempty()          { return m_symmX.empty(); }
  bool symmYempty()          { return m_symmY.empty(); }
  bool symmZempty()          { return m_symmZ.empty(); }

  Real_t& arealg(Index_t idx)     { return m_arealg[idx] ; }
  Real_t& ss(Index_t idx)         { return m_ss[idx] ; }

  Index_t nodeElemCount(Index_t idx)
  { return m_nodeElemStart[idx+1] - m_nodeElemStart[idx] ; }

  Index_t *nodeElemCornerList(Index_t idx)
  { return &m_nodeElemCornerList[m_nodeElemStart[idx]] ; }

  Index_t numElem() const           { return m_nElem[0]*m_nElem[1]*m_nElem[2]; }
  Index_t numElem(size_t dim) const { return m_nElem[dim]; }

  Index_t numNode() const           { return m_nNode[0]*m_nNode[1]*m_nNode[2]; }
  Index_t numNode(size_t dim) const { return m_nNode[dim]; }

  // for compatibility with the old code
  Index_t sizeX() const          { return m_nElem[0]; }
  Index_t sizeY() const          { return m_nElem[1]; }
  Index_t sizeZ() const          { return m_nElem[2]; }

  Index_t colLoc() const         { return m_ts.z(dash::myid()); }
  Index_t rowLoc() const         { return m_ts.y(dash::myid()); }
  Index_t planeLoc() const       { return m_ts.x(dash::myid()); }
  Index_t tp(size_t dim) const   { return m_ts.extent(dim);     }
  Index_t tp() const             { return m_ts.extent(0);       }

  Index_t numRanks() const       { return dash::size(); }
  Index_t maxEdgeSize() const    { return 1+std::max({sizeX(), sizeY(), sizeZ()}); }
  Index_t maxPlaneSize() const   { return maxEdgeSize()*maxEdgeSize(); }


  // Cutoffs
  Real_t u_cut() const               { return m_param.m_u_cut; }
  Real_t e_cut() const               { return m_param.m_e_cut; }
  Real_t p_cut() const               { return m_param.m_p_cut; }
  Real_t q_cut() const               { return m_param.m_q_cut; }
  Real_t v_cut() const               { return m_param.m_v_cut; }

  // Other constants (usually are settable via input file in real codes)
  Real_t hgcoef() const              { return m_param.m_hgcoef; }
  Real_t qstop() const               { return m_param.m_qstop; }
  Real_t monoq_max_slope() const     { return m_param.m_monoq_max_slope; }
  Real_t monoq_limiter_mult() const  { return m_param.m_monoq_limiter_mult; }
  Real_t ss4o3() const               { return m_param.m_ss4o3; }
  Real_t qlc_monoq() const           { return m_param.m_qlc_monoq; }
  Real_t qqc_monoq() const           { return m_param.m_qqc_monoq; }
  Real_t qqc() const                 { return m_param.m_qqc; }

  Real_t eosvmax() const             { return m_param.m_eosvmax; }
  Real_t eosvmin() const             { return m_param.m_eosvmin; }
  Real_t pmin() const                { return m_param.m_pmin; }
  Real_t emin() const                { return m_param.m_emin; }
  Real_t dvovmax() const             { return m_param.m_dvovmax; }
  Real_t refdens() const             { return m_param.m_refdens; }

  // Timestep controls, etc...
  Real_t& time()                     { return m_param.m_time; }
  Real_t& deltatime()                { return m_param.m_deltatime; }
  Real_t& deltatimemultlb()          { return m_param.m_deltatimemultlb; }
  Real_t& deltatimemultub()          { return m_param.m_deltatimemultub; }
  Real_t& stoptime()                 { return m_param.m_stoptime; }
  Real_t& dtcourant()                { return m_param.m_dtcourant; }
  Real_t& dthydro()                  { return m_param.m_dthydro; }
  Real_t& dtmax()                    { return m_param.m_dtmax; }
  Real_t& dtfixed()                  { return m_param.m_dtfixed; }
  Int_t&  cycle()                    { return m_param.m_cycle; }

  Int_t&  chunksize()                { return m_chunksize; }


  Index_t&  numReg()                 { return m_region.numReg(); }
  Int_t&    cost()                   { return m_region.cost(); }
  Index_t&  regElemSize(Index_t idx) { return m_region.regElemSize(idx); }
  Index_t&  regNumList(Index_t idx)  { return m_region.regNumList(idx); }
  Index_t*  regNumList()             { return m_region.regNumList(); }
  Index_t*  regElemlist(Int_t r)     { return m_region.regElemlist(r); }
  Index_t&  regElemlist(Int_t r, Index_t idx) { return m_region.regElemlist(r,idx); }

  void InitialBoundaryExchange();

  void TimeIncrement();
  void LagrangeLeapFrog();
  void LagrangeNodal();
  void LagrangeElements();

  void CalcForceForNodes();
  void CalcAccelerationForNodes(Index_t numNode);
  void ApplyAccelerationBoundaryConditionsForNodes();

  void CalcVelocityForNodes(const Real_t dt,
			    const Real_t u_cut, Index_t numNode);
  void CalcPositionForNodes(const Real_t dt, Index_t numNode);

  void CalcLagrangeElements();
  void CalcQForElems();
  void ApplyMaterialPropertiesForElems();
  void UpdateVolumesForElems(Real_t v_cut, Index_t length);

  void VerifyAndWriteFinalOutput(Real_t elapsed,
				 Int_t  nx,
				 Int_t  numRanks);
};



#endif /* LULESH_DASH_H_INCLUDED */

