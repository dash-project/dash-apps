
#include <iostream>
#include "lulesh-util.h"

void print_config(Domain& d,
		  std::ostream& os)
{
  auto myid = dash::myid();
  os << "[ " << myid << "] colLoc()             : " << d.colLoc()       << std::endl;
  os << "[ " << myid << "] rowLoc()             : " << d.rowLoc()       << std::endl;
  os << "[ " << myid << "] planeLoc()           : " << d.planeLoc()     << std::endl;
  os << "[ " << myid << "] numElem()            : " << d.numElem()      << std::endl;
  os << "[ " << myid << "] numElem(0)           : " << d.numElem(0)     << std::endl;
  os << "[ " << myid << "] numElem(1)           : " << d.numElem(1)     << std::endl;
  os << "[ " << myid << "] numElem(2)           : " << d.numElem(2)     << std::endl;
  os << "[ " << myid << "] numNode()            : " << d.numNode()      << std::endl;
  os << "[ " << myid << "] numNode(0)           : " << d.numNode(0)     << std::endl;
  os << "[ " << myid << "] numNode(1)           : " << d.numNode(1)     << std::endl;
  os << "[ " << myid << "] numNode(2)           : " << d.numNode(2)     << std::endl;
  os << "[ " << myid << "] sizeX()              : " << d.sizeX()        << std::endl;
  os << "[ " << myid << "] sizeY()              : " << d.sizeY()        << std::endl;
  os << "[ " << myid << "] sizeZ()              : " << d.sizeZ()        << std::endl;
  os << "[ " << myid << "] maxEdgeSize()        : " << d.maxEdgeSize()  << std::endl;
  os << "[ " << myid << "] maxPlaneSize()       : " << d.maxPlaneSize() << std::endl;

  os << "[ " << myid << "] u_cut()              : " << d.u_cut() << std::endl;
  os << "[ " << myid << "] e_cut()              : " << d.e_cut() << std::endl;
  os << "[ " << myid << "] p_cut()              : " << d.p_cut() << std::endl;
  os << "[ " << myid << "] q_cut()              : " << d.q_cut() << std::endl;
  os << "[ " << myid << "] v_cut()              : " << d.v_cut() << std::endl;

  os << "[ " << myid << "] hgcoef()             : " << d.hgcoef()             << std::endl;
  os << "[ " << myid << "] qstop()              : " << d.qstop()              << std::endl;
  os << "[ " << myid << "] monoq_max_slope()    : " << d.monoq_max_slope()    << std::endl;
  os << "[ " << myid << "] monoq_limiter_mult() : " << d.monoq_limiter_mult() << std::endl;
  os << "[ " << myid << "] ss4o3() const        : " << d.ss4o3()              << std::endl;
  os << "[ " << myid << "] qlc_monoq() const    : " << d.qlc_monoq()          << std::endl;
  os << "[ " << myid << "] qqc_monoq() const    : " << d.qqc_monoq()          << std::endl;
  os << "[ " << myid << "] qqc() const          : " << d.qqc()                << std::endl;

  os << "[ " << myid << "] eosvmax() const      : " << d.eosvmax() << std::endl;
  os << "[ " << myid << "] eosvmin() const      : " << d.eosvmin() << std::endl;
  os << "[ " << myid << "] pmin() const         : " << d.pmin()    << std::endl;
  os << "[ " << myid << "] emin() const         : " << d.emin()    << std::endl;
  os << "[ " << myid << "] dvovmax() const      : " << d.dvovmax() << std::endl;
  os << "[ " << myid << "] refdens() const      : " << d.refdens() << std::endl;

  os << "[ " << myid << "] time()               : " << d.time()            << std::endl;
  os << "[ " << myid << "] deltatime()          : " << d.deltatime()       << std::endl;
  os << "[ " << myid << "] deltatimemultlb()    : " << d.deltatimemultlb() << std::endl;
  os << "[ " << myid << "] deltatimemultub()    : " << d.deltatimemultub() << std::endl;
  os << "[ " << myid << "] stoptime()           : " << d.stoptime()        << std::endl;
  os << "[ " << myid << "] dtcourant()          : " << d.dtcourant()       << std::endl;
  os << "[ " << myid << "] dthydro()            : " << d.dthydro()         << std::endl;
  os << "[ " << myid << "] dtmax()              : " << d.dtmax()           << std::endl;
  os << "[ " << myid << "] dtfixed()            : " << d.dtfixed()         << std::endl;
  os << "[ " << myid << "] cycle()              : " << d.cycle()           << std::endl;
}
