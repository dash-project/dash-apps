
#include <iostream>

#include "lulesh-dash.h"
#include "lulesh-comm-dash.h"
#include "lulesh-comm-dash-onesided.h"

#ifdef DASH_USE_GET

using namespace std;

char tagname[][8] =
  {
    "X0     ",
    "X1     ",
    "Y0     ",
    "Y1     ",
    "Z0     ",
    "Z1     ",
    "X0Y0   ",
    "X0Y1   ",
    "X1Y0   ",
    "X1Y1   ",
    "X0Z0   ",
    "X0Z1   ",
    "X1Z0   ",
    "X1Z1   ",
    "Y0Z0   ",
    "Y0Z1   ",
    "Y1Z0   ",
    "Y1Z1   ",
    "X0Y0Z0 ",
    "X0Y0Z1 ",
    "X0Y1Z0 ",
    "X0Y1Z1 ",
    "X1Y0Z0 ",
    "X1Y0Z1 ",
    "X1Y1Z0 ",
    "X1Y1Z1 ",
  };


static void
get_yield(const dash::GlobIter<double, dash::BlockPattern<1> > src,
          Real_t *srcAddr, size_t recvCount)
{
  auto fut =
    dash::copy_async(
      src,
      src + recvCount,
      srcAddr);
  while(!fut.test()) dash::tasks::yield();
}

// debug output for syncs performed below
void DBGSYNC(int fields, int elem, int tag)
{
  /*
  cout << "DBGSYNC [" << dash::myid() << "] " << elem
       << " " << tagname[tag] << endl;
  */
}

void DASHCommPut(Domain& domain, DASHComm& comm,
                 Index_t xferFields, Domain_member *fieldData,
                 Index_t dx, Index_t dy, Index_t dz,
                 bool doSend, bool planeOnly)
{
  Real_t *destAddr;
  bool rowNotMin, rowNotMax, colNotMin, colNotMax, planeNotMin, planeNotMax;

  Index_t maxPlaneComm = xferFields * domain.maxPlaneSize();
  Index_t maxEdgeComm  = xferFields * domain.maxEdgeSize();

  Index_t pmsg = 0; // plane comm msg
  Index_t emsg = 0; // edge comm msg
  Index_t cmsg = 0; // corner comm msg

  // assume communication to 6 neighbors by default
  rowNotMin = rowNotMax = colNotMin = colNotMax = planeNotMin = planeNotMax = true;

  if( domain.rowLoc()   == 0 )               { rowNotMin   = false; }
  if( domain.rowLoc()   == (domain.tp()-1) ) { rowNotMax   = false; }
  if( domain.colLoc()   == 0 )               { colNotMin   = false; }
  if( domain.colLoc()   == (domain.tp()-1) ) { colNotMax   = false; }
  if( domain.planeLoc() == 0 )               { planeNotMin = false; }
  if( domain.planeLoc() == (domain.tp()-1) ) { planeNotMax = false; }

  int myRank = dash::myid();

  if( planeNotMin | planeNotMax ) {
    // ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE

    if (planeNotMin) {
      auto dest = &comm.commDataSend()[comm.offset(Z0, xferFields)];
      dash::tasks::ASYNC(
        [=, &domain](){
          int sendCount = dx * dy;
          Real_t *destAddr = dest;
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<sendCount; ++i) {
              destAddr[i] = (domain.*src)(i);
            }
            destAddr += sendCount;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );
      /*
        MPI_Isend(destAddr, xferFields*sendCount, baseType,
        myRank - domain.tp()*domain.tp(), msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg]);
      */

      ++pmsg;
    }
    if (planeNotMax && doSend) {
      //auto dest = comm.dest(myRank + (domain.tp() * domain.tp()), Z0);
      auto dest = &comm.commDataSend()[comm.offset(Z1, xferFields)];
      dash::tasks::ASYNC(
        [=, &domain](){
          int sendCount = dx * dy;
          Real_t *destAddr = dest;
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<sendCount; ++i) {
              destAddr[i] = (domain.*src)(dx*dy*(dz - 1) + i);
            }
            destAddr += sendCount;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );
      /*
        MPI_Isend(destAddr, xferFields*sendCount, baseType,
        myRank + domain.tp()*domain.tp(), msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg]);
      */
      ++pmsg;
    }
  }
  if (rowNotMin | rowNotMax) {
    // ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE
    if (rowNotMin) {
      auto dest = &comm.commDataSend()[comm.offset(Y0, xferFields)];
      dash::tasks::ASYNC(
        [=, &domain](){
          int sendCount = dx * dz;
          Real_t *destAddr = dest;
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              for (Index_t j=0; j<dx; ++j) {
                destAddr[i*dx+j] = (domain.*src)(i*dx*dy + j);
              }
            }
            destAddr += sendCount;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );
      /*
        MPI_Isend(destAddr, xferFields*sendCount, baseType,
        myRank - domain.tp(), msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg]);
      */
      ++pmsg;
    }
    if (rowNotMax && doSend) {
      auto dest = &comm.commDataSend()[comm.offset(Y1, xferFields)];
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          int sendCount = dx * dz;
          Real_t *destAddr = dest;
          destAddr = &comm.commDataSend()[pmsg * maxPlaneComm];
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              for (Index_t j=0; j<dx; ++j) {
                destAddr[i*dx+j] = (domain.*src)(dx*(dy - 1) + i*dx*dy + j);
              }
            }
            destAddr += sendCount;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );

      /*
        MPI_Isend(destAddr, xferFields*sendCount, baseType,
        myRank + domain.tp(), msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg]);
      */
      ++pmsg;
    }
  }
  if (colNotMin | colNotMax) {
    // ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE

    if (colNotMin) {
      auto dest = &comm.commDataSend()[comm.offset(X0, xferFields)];
      dash::tasks::ASYNC(
        [=, &domain](){
          int sendCount = dy * dz;
          Real_t *destAddr = dest;
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              for (Index_t j=0; j<dy; ++j) {
                destAddr[i*dy + j] = (domain.*src)(i*dx*dy + j*dx);
              }
            }
            destAddr += sendCount;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );

      /*
        MPI_Isend(destAddr, xferFields*sendCount, baseType,
        myRank - 1, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg]);
      */
      ++pmsg;
    }
    if (colNotMax && doSend) {
      auto dest = &comm.commDataSend()[comm.offset(X1, xferFields)];
      dash::tasks::ASYNC(
        [=, &domain](){
          int sendCount = dy * dz;
          Real_t *destAddr = dest;
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              for (Index_t j=0; j<dy; ++j) {
                destAddr[i*dy + j] = (domain.*src)(dx - 1 + i*dx*dy + j*dx);
              }
            }
            destAddr += sendCount;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );

      /*
        MPI_Isend(destAddr, xferFields*sendCount, baseType,
        myRank + 1, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg]);
      */
      ++pmsg;
    }
  }

  if (!planeOnly) {
    if (rowNotMin && colNotMin) {
      auto dest  = &comm.commDataSend()[comm.offset(X0Y0, xferFields)];
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          Real_t *destAddr = dest;
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              destAddr[i] = (domain.*src)(i*dx*dy);
            }
            destAddr += dz;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );

      /*
        MPI_Isend(destAddr, xferFields*dz, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (rowNotMin && planeNotMin) {
      auto dest  = &comm.commDataSend()[comm.offset(Y0Z0, xferFields)];
      dash::tasks::ASYNC(
        [=, &domain](){
          Real_t *destAddr = dest;
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dx; ++i) {
              destAddr[i] = (domain.*src)(i);
            }
            destAddr += dx;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );

      /*
        MPI_Isend(destAddr, xferFields*dx, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (colNotMin && planeNotMin) {
      auto dest  = &comm.commDataSend()[comm.offset(X0Z0, xferFields)];
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          Real_t *destAddr = dest;
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dy; ++i) {
              destAddr[i] = (domain.*src)(i*dx);
            }
            destAddr += dy;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );

      /*
        MPI_Isend(destAddr, xferFields*dy, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (rowNotMax && colNotMax && doSend) {
      auto dest  = &comm.commDataSend()[comm.offset(X1Y1, xferFields)];
      dash::tasks::ASYNC(
        [=, &domain](){
          Real_t *destAddr = dest;
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              destAddr[i] = (domain.*src)(dx*dy - 1 + i*dx*dy);
            }
            destAddr += dz;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );
      /*
        MPI_Isend(destAddr, xferFields*dz, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (rowNotMax && planeNotMax && doSend) {
      auto dest  = &comm.commDataSend()[comm.offset(Y1Z1, xferFields)];
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          Real_t *destAddr = dest;
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dx; ++i) {
              destAddr[i] = (domain.*src)(dx*(dy-1) + dx*dy*(dz-1) + i);
            }
            destAddr += dx;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );
      /*
        MPI_Isend(destAddr, xferFields*dx, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (colNotMax && planeNotMax && doSend) {
      auto dest  = &comm.commDataSend()[comm.offset(X1Z1, xferFields)];
      dash::tasks::ASYNC(
        [=, &domain](){
          Real_t *destAddr = dest;
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dy; ++i) {
              destAddr[i] = (domain.*src)(dx*dy*(dz-1) + dx - 1 + i*dx);
            }
            destAddr += dy;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );
      /*
        MPI_Isend(destAddr, xferFields*dy, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (rowNotMax && colNotMin && doSend) {
      auto dest  = &comm.commDataSend()[comm.offset(X0Y1, xferFields)];
      dash::tasks::ASYNC(
        [=, &domain](){
          Real_t *destAddr = dest;
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              destAddr[i] = (domain.*src)(dx*(dy-1) + i*dx*dy);
            }
            destAddr += dz;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );
      /*
        MPI_Isend(destAddr, xferFields*dz, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (rowNotMin && planeNotMax && doSend) {
      auto dest  = &comm.commDataSend()[comm.offset(Y0Z1, xferFields)];
      dash::tasks::ASYNC(
        [=, &domain](){
          Real_t *destAddr = dest;
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dx; ++i) {
              destAddr[i] = (domain.*src)(dx*dy*(dz-1) + i);
            }
            destAddr += dx;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );
      /*
        MPI_Isend(destAddr, xferFields*dx, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (colNotMin && planeNotMax && doSend) {
      auto dest  = &comm.commDataSend()[comm.offset(X0Z1, xferFields)];
      dash::tasks::ASYNC(
        [=, &domain](){
          Real_t *destAddr = dest;
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dy; ++i) {
              destAddr[i] = (domain.*src)(dx*dy*(dz-1) + i*dx);
            }
            destAddr += dy;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );
      /*
        MPI_Isend(destAddr, xferFields*dy, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (rowNotMin && colNotMax) {
      auto dest  = &comm.commDataSend()[comm.offset(X1Y0, xferFields)];
      dash::tasks::ASYNC(
        [=, &domain](){
          Real_t *destAddr = dest;
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              destAddr[i] = (domain.*src)(dx - 1 + i*dx*dy);
            }
            destAddr += dz;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );
      /*
        MPI_Isend(destAddr, xferFields*dz, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (rowNotMax && planeNotMin) {
      auto dest  = &comm.commDataSend()[comm.offset(Y1Z0, xferFields)];
      dash::tasks::ASYNC(
        [=, &domain](){
          Real_t *destAddr = dest;
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dx; ++i) {
              destAddr[i] = (domain.*src)(dx*(dy - 1) + i);
            }
            destAddr += dx;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );

      /*
        MPI_Isend(destAddr, xferFields*dx, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (colNotMax && planeNotMin) {
      auto dest  = &comm.commDataSend()[comm.offset(X1Z0, xferFields)];
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          Real_t *destAddr = dest;
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dy; ++i) {
              destAddr[i] = (domain.*src)(dx - 1 + i*dx);
            }
            destAddr += dy;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );

      /*
        MPI_Isend(destAddr, xferFields*dy, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (rowNotMin && colNotMin && planeNotMin) {
      auto dest = &comm.commDataSend()[X0Y0Z0];
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          // corner at domain logical coord (0, 0, 0)
          Real_t *comBuf = dest;
          for (Index_t fi=0; fi<xferFields; ++fi) {
            comBuf[fi] = (domain.*fieldData[fi])(0);
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );

      /*
        MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]);
      */
      ++cmsg;
    }
    if (rowNotMin && colNotMin && planeNotMax && doSend) {
      auto dest = &comm.commDataSend()[comm.offset(X0Y0Z1, xferFields)];
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          // corner at domain logical coord (0, 0, 1)
          Real_t *comBuf = dest;
          Index_t idx = dx*dy*(dz - 1);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            comBuf[fi] = (domain.*fieldData[fi])(idx);
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );

      /*
        MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]);
      */
      ++cmsg;
    }
    if (rowNotMin && colNotMax && planeNotMin) {
      auto dest = &comm.commDataSend()[comm.offset(X1Y0Z0, xferFields)];
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          // corner at domain logical coord (1, 0, 0)
          Real_t *comBuf = dest;
          Index_t idx = dx - 1;
          for (Index_t fi=0; fi<xferFields; ++fi) {
            comBuf[fi] = (domain.*fieldData[fi])(idx);
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );

      /*
        MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]);
      */
      ++cmsg;
    }
    if (rowNotMin && colNotMax && planeNotMax && doSend) {
      auto dest = &comm.commDataSend()[X1Y0Z1];
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          // corner at domain logical coord (1, 0, 1)
          Real_t *comBuf = dest;
          Index_t idx = dx*dy*(dz - 1) + (dx - 1);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            comBuf[fi] = (domain.*fieldData[fi])(idx);
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );

      /*
        MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]);
      */
      ++cmsg;
    }
    if (rowNotMax && colNotMin && planeNotMin) {
      auto dest = &comm.commDataSend()[comm.offset(X0Y1Z0, xferFields)];
      dash::tasks::ASYNC(
        [=, &domain](){
          // corner at domain logical coord (0, 1, 0)
          Real_t *comBuf = dest;
          Index_t idx = dx*(dy - 1);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            comBuf[fi] = (domain.*fieldData[fi])(idx);
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );

      /*
        MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]);
      */
      ++cmsg;
    }
    if (rowNotMax && colNotMin && planeNotMax && doSend) {
      auto dest = &comm.commDataSend()[comm.offset(X0Y1Z1, xferFields)];
      dash::tasks::ASYNC(
        [=, &domain](){
          // corner at domain logical coord (0, 1, 1)
          Real_t *comBuf = dest;
          Index_t idx = dx*dy*(dz - 1) + dx*(dy - 1);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            comBuf[fi] = (domain.*fieldData[fi])(idx);
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );

      /*
        MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]);
      */
      ++cmsg;
    }
    if (rowNotMax && colNotMax && planeNotMin) {
      auto dest = &comm.commDataSend()[comm.offset(X1Y1Z0, xferFields)];
      dash::tasks::ASYNC(
        [=, &domain](){
          // corner at domain logical coord (1, 1, 0)
          Real_t *comBuf = dest;
          Index_t idx = dx*dy - 1;
          for (Index_t fi=0; fi<xferFields; ++fi) {
            comBuf[fi] = (domain.*fieldData[fi])(idx);
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );

      /*
        MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]);
      */
      ++cmsg;
    }
    if (rowNotMax && colNotMax && planeNotMax && doSend) {
      auto dest = &comm.commDataSend()[comm.offset(X1Y1Z1, xferFields)];
      dash::tasks::ASYNC(
        [=, &domain](){
          // corner at domain logical coord (1, 1, 1)
          Real_t *comBuf = dest;
          Index_t idx = dx*dy*dz - 1;
          for (Index_t fi=0; fi<xferFields; ++fi) {
            comBuf[fi] = (domain.*fieldData[fi])(idx);
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::out(dest)
      );

      /*
        MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]);
      */
      ++cmsg;
    }
  }
}

void DASHCommSBN(Domain& domain, DASHComm& comm, int xferFields,
                 Domain_member *fieldData)
{
  if (domain.numRanks() == 1)
    return;

  /* summation order should be from smallest value to largest */
  /* or we could try out kahan summation! */

  int myRank = dash::myid();

  Index_t dx = domain.sizeX() + 1;
  Index_t dy = domain.sizeY() + 1;
  Index_t dz = domain.sizeZ() + 1;

  Index_t rowNotMin, rowNotMax, colNotMin, colNotMax, planeNotMin, planeNotMax;

  // assume communication to 6 neighbors by default
  rowNotMin = rowNotMax = colNotMin = colNotMax = planeNotMin = planeNotMax = true;

  if( domain.rowLoc()   == 0 )               { rowNotMin   = false; }
  if( domain.rowLoc()   == (domain.tp()-1) ) { rowNotMax   = false; }
  if( domain.colLoc()   == 0 )               { colNotMin   = false; }
  if( domain.colLoc()   == (domain.tp()-1) ) { colNotMax   = false; }
  if( domain.planeLoc() == 0 )               { planeNotMin = false; }
  if( domain.planeLoc() == (domain.tp()-1) ) { planeNotMax = false; }

  if (planeNotMin | planeNotMax) {
    /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
    Index_t opCount = dx * dy;

    if (planeNotMin) {
      auto src = comm.src(myRank - (domain.tp() * domain.tp()), Z1, xferFields);
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(Z1, xferFields)];
          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          DBGSYNC(xferFields, opCount, Z1);

          size_t recvCount = xferFields*opCount;
          get_yield(src, srcAddr, recvCount);

          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<opCount; ++i) {
              (domain.*dest)(i) += srcAddr[i];
            }
            srcAddr += opCount;
          }

        },
        dash::tasks::in(src),
        // TODO: dummy IN dependencies, replace them with CONCURRENT!!
        dash::tasks::in(&(domain.*fieldData[0])(0))
      );
    }
    if (planeNotMax) {
      auto src = comm.src(myRank + (domain.tp() * domain.tp()), Z0, xferFields);
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(Z0, xferFields)];
          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          DBGSYNC(xferFields, opCount, Z0);

          size_t recvCount = xferFields*opCount;
          get_yield(src, srcAddr, recvCount);

          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<opCount; ++i) {
              (domain.*dest)(dx*dy*(dz - 1) + i) += srcAddr[i];
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(src),
        // TODO: dummy IN dependencies, replace them with CONCURRENT!!
        dash::tasks::in(&(domain.*fieldData[0])(0))
      );
    }
  }

  if (rowNotMin | rowNotMax) {
    /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
    Index_t opCount = dx * dz;

    if (rowNotMin) {
      auto src = comm.src(myRank - domain.tp(), Y1, xferFields);
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(Y1, xferFields)];
          DBGSYNC(xferFields, opCount, Y1);

          size_t recvCount = xferFields*opCount;
          get_yield(src, srcAddr, recvCount);

          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              for (Index_t j=0; j<dx; ++j) {
                (domain.*dest)(i*dx*dy + j) += srcAddr[i*dx + j];
              }
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(src),
        // TODO: dummy IN dependencies, replace them with CONCURRENT!!
        dash::tasks::in(&(domain.*fieldData[0])(0))
      );
    }
    if (rowNotMax) {
      auto src = comm.src(myRank + domain.tp(), Y0, xferFields);
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(Y0, xferFields)];
          DBGSYNC(xferFields, opCount, Y0);

          size_t recvCount = xferFields*opCount;
          get_yield(src, srcAddr, recvCount);

          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              for (Index_t j=0; j<dx; ++j) {
                (domain.*dest)(dx*(dy - 1) + i*dx*dy + j) += srcAddr[i*dx + j];
              }
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(src),
        // TODO: dummy IN dependencies, replace them with CONCURRENT!!
        dash::tasks::in(&(domain.*fieldData[0])(0))
      );
    }
  }
  if (colNotMin | colNotMax) {
    /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
    Index_t opCount = dy * dz;

    if (colNotMin) {
      auto src = comm.src( myRank - 1, X1, xferFields);
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(X1, xferFields)];
          DBGSYNC(xferFields, opCount, X1);

          size_t recvCount = xferFields*opCount;
          get_yield(src, srcAddr, recvCount);

          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              for (Index_t j=0; j<dy; ++j) {
                (domain.*dest)(i*dx*dy + j*dx) += srcAddr[i*dy + j];
              }
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(src),
        // TODO: dummy IN dependencies, replace them with CONCURRENT!!
        dash::tasks::in(&(domain.*fieldData[0])(0))
      );
    }
    if (colNotMax) {
      auto src = comm.src( myRank + 1, X0, xferFields );
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(X0, xferFields)];
          DBGSYNC(xferFields, opCount, X0);

          size_t recvCount = xferFields*opCount;
          get_yield(src, srcAddr, recvCount);

          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              for (Index_t j=0; j<dy; ++j) {
                (domain.*dest)(dx - 1 + i*dx*dy + j*dx) += srcAddr[i*dy + j];
              }
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(src),
        // TODO: dummy IN dependencies, replace them with CONCURRENT!!
        dash::tasks::in(&(domain.*fieldData[0])(0))
      );
    }
  }

  if (rowNotMin & colNotMin) {
    int rank = myRank - domain.tp() - 1;
    auto src  = comm.src( rank, X1Y1, xferFields);
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X1Y1, xferFields)];
        DBGSYNC(xferFields, dz, X1Y1);

        size_t recvCount = xferFields*dz;
        get_yield(src, srcAddr, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dz; ++i) {
            (domain.*dest)(i*dx*dy) += srcAddr[i];
          }
          srcAddr += dz;
        }
      },
      dash::tasks::in(src),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }

  if (rowNotMin & planeNotMin) {
    int  rank = myRank - domain.tp()*domain.tp() - domain.tp();
    auto src  = comm.src( rank, Y1Z1, xferFields);
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(Y1Z1, xferFields)];
        DBGSYNC(xferFields, dx, Y1Z1);

        size_t recvCount = xferFields*dx;
        get_yield(src, srcAddr, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dx; ++i) {
            (domain.*dest)(i) += srcAddr[i];
          }
          srcAddr += dx;
        }
      },
      dash::tasks::in(src),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }

  if (colNotMin & planeNotMin) {
    int  rank = myRank - domain.tp()*domain.tp() - 1;
    auto src  = comm.src( rank, X1Z1, xferFields);
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X1Z1, xferFields)];
        DBGSYNC(xferFields, dy, X1Z1);

        size_t recvCount = xferFields*dy;
        get_yield(src, srcAddr, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dy; ++i) {
            (domain.*dest)(i*dx) += srcAddr[i];
          }
          srcAddr += dy;
        }
      },
      dash::tasks::in(src),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }

  if (rowNotMax & colNotMax) {
    int  rank = myRank + domain.tp() + 1;
    auto src  = comm.src( rank, X0Y0, xferFields);
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X0Y0, xferFields)];
        DBGSYNC(xferFields, dz, X0Y0);

        size_t recvCount = xferFields*dz;
        get_yield(src, srcAddr, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dz; ++i) {
            (domain.*dest)(dx*dy - 1 + i*dx*dy) += srcAddr[i];
          }
          srcAddr += dz;
        }
      },
      dash::tasks::in(src),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }

  if (rowNotMax & planeNotMax) {
    int  rank = myRank + domain.tp()*domain.tp() + domain.tp();
    auto src  = comm.src( rank, Y0Z0, xferFields);
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(Y0Z0, xferFields)];
        DBGSYNC(xferFields, dx, Y0Z0);

        size_t recvCount = xferFields*dx;
        get_yield(src, srcAddr, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dx; ++i) {
            (domain.*dest)(dx*(dy-1) + dx*dy*(dz-1) + i) += srcAddr[i];
          }
          srcAddr += dx;
        }
      },
      dash::tasks::in(src),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }

  if (colNotMax & planeNotMax) {
    int  rank = myRank + domain.tp()*domain.tp() + 1;
    auto src  = comm.src( rank, X0Z0, xferFields);
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X0Z0, xferFields)];
        DBGSYNC(xferFields, dy, X0Z0);

        size_t recvCount = xferFields*dy;
        get_yield(src, srcAddr, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dy; ++i) {
            (domain.*dest)(dx*dy*(dz-1) + dx - 1 + i*dx) += srcAddr[i];
          }
          srcAddr += dy;
        }
      },
      dash::tasks::in(src),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }

  if (rowNotMax & colNotMin) {
    int  rank = myRank + domain.tp() - 1;
    auto src  = comm.src( rank, X1Y0, xferFields);
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X1Y0, xferFields)];
        DBGSYNC(xferFields, dz, X1Y0);

        size_t recvCount = xferFields*dz;
        get_yield(src, srcAddr, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dz; ++i) {
            (domain.*dest)(dx*(dy-1) + i*dx*dy) += srcAddr[i];
          }
          srcAddr += dz;
        }
      },
      dash::tasks::in(src),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }

  if (rowNotMin & planeNotMax) {
    int  rank = myRank + domain.tp()*domain.tp() - domain.tp();
    auto src  = comm.src( rank, Y1Z0, xferFields);
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(Y1Z0, xferFields)];
        DBGSYNC(xferFields, dx, Y1Z0);

        size_t recvCount = xferFields*dx;
        get_yield(src, srcAddr, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dx; ++i) {
            (domain.*dest)(dx*dy*(dz-1) + i) += srcAddr[i];
          }
          srcAddr += dx;
        }
      },
      dash::tasks::in(src),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }

  if (colNotMin & planeNotMax) {
    int  rank = myRank + domain.tp()*domain.tp() - 1;
    auto src  = comm.src( rank, X1Z0, xferFields);
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X1Z0, xferFields)];
        DBGSYNC(xferFields, dy, X1Z0);

        size_t recvCount = xferFields*dy;
        get_yield(src, srcAddr, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dy; ++i) {
            (domain.*dest)(dx*dy*(dz-1) + i*dx) += srcAddr[i];
          }
          srcAddr += dy;
        }
      },
      dash::tasks::in(src),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }

  if (rowNotMin & colNotMax) {
    int  rank = myRank - domain.tp() + 1;
    auto src  = comm.src( rank, X0Y1, xferFields);
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X0Y1, xferFields)];
        DBGSYNC(xferFields, dz, X0Y1);

        size_t recvCount = xferFields*dz;
        get_yield(src, srcAddr, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dz; ++i) {
            (domain.*dest)(dx - 1 + i*dx*dy) += srcAddr[i];
          }
          srcAddr += dz;
        }
      },
      dash::tasks::in(src),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }

  if (rowNotMax & planeNotMin) {
    int  rank = myRank - domain.tp()*domain.tp() + domain.tp();
    auto src  = comm.src( rank, Y0Z1, xferFields );
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(Y0Z1, xferFields)];
        DBGSYNC(xferFields, dx, Y0Z1);

        size_t recvCount = xferFields*dx;
        get_yield(src, srcAddr, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dx; ++i) {
            (domain.*dest)(dx*(dy - 1) + i) += srcAddr[i];
          }
          srcAddr += dx;
        }
      },
      dash::tasks::in(src),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }

  if (colNotMax & planeNotMin) {
    int  rank = myRank - domain.tp()*domain.tp() + 1;
    auto src  = comm.src( rank, X0Z1, xferFields );
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X0Z1, xferFields)];
        DBGSYNC(xferFields, dy, X0Z1);

        size_t recvCount = xferFields*dy;
        get_yield(src, srcAddr, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dy; ++i) {
            (domain.*dest)(dx - 1 + i*dx) += srcAddr[i];
          }
          srcAddr += dy;
        }
      },
      dash::tasks::in(src),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }

  if (rowNotMin & colNotMin & planeNotMin) {
    int  rank = myRank - domain.tp()*domain.tp() - domain.tp() - 1;
    auto src  = comm.dest( rank, X1Y1Z1, xferFields );
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        /* corner at domain logical coord (0, 0, 0) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X1Y1Z1, xferFields)];
        DBGSYNC(xferFields, 1, X1Y1Z1);

        size_t recvCount = xferFields;
        get_yield(src, comBuf, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(0) += comBuf[fi];
        }
      },
      dash::tasks::in(src),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }
  if (rowNotMin & colNotMin & planeNotMax) {
    int  rank = myRank + domain.tp()*domain.tp() - domain.tp() - 1;
    auto src  = comm.src( rank, X1Y1Z0, xferFields );
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        /* corner at domain logical coord (0, 0, 1) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X1Y1Z0, xferFields)];
        DBGSYNC(xferFields, 1, X1Y1Z0);
        Index_t idx = dx*dy*(dz - 1);

        size_t recvCount = xferFields;
        get_yield(src, comBuf, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) += comBuf[fi];
        }
      },
      dash::tasks::in(src),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }
  if (rowNotMin & colNotMax & planeNotMin) {
    int  rank = myRank - domain.tp()*domain.tp() - domain.tp() + 1;
    auto src  = comm.src( rank, X0Y1Z1, xferFields );
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        /* corner at domain logical coord (1, 0, 0) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X0Y1Z1, xferFields)];
        DBGSYNC(xferFields, 1, X0Y1Z1);
        Index_t idx = dx - 1;

        size_t recvCount = xferFields;
        get_yield(src, comBuf, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) += comBuf[fi];
        }
      },
      dash::tasks::in(src),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }
  if (rowNotMin & colNotMax & planeNotMax) {
    int  rank = myRank + domain.tp()*domain.tp() - domain.tp() + 1;
    auto src  = comm.dest( rank, X0Y1Z0, xferFields );
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        /* corner at domain logical coord (1, 0, 1) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X0Y1Z0, xferFields)];
        DBGSYNC(xferFields, 1, X0Y1Z0);
        Index_t idx = dx*dy*(dz - 1) + (dx - 1);

        size_t recvCount = xferFields;
        get_yield(src, comBuf, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) += comBuf[fi];
        }
      },
      dash::tasks::in(src),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }
  if (rowNotMax & colNotMin & planeNotMin) {
    int  rank = myRank - domain.tp()*domain.tp() + domain.tp() - 1;
    auto src  = comm.src( rank, X1Y0Z1, xferFields );
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        /* corner at domain logical coord (0, 1, 0) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X1Y0Z1, xferFields)];
        DBGSYNC(xferFields, 1, X1Y0Z1);
        Index_t idx = dx*(dy - 1);

        size_t recvCount = xferFields;
        get_yield(src, comBuf, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) += comBuf[fi];
        }
      },
      dash::tasks::in(src),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }
  if (rowNotMax & colNotMin & planeNotMax) {
    int  rank = myRank + domain.tp()*domain.tp() + domain.tp() - 1;
    auto src  = comm.src( rank, X1Y0Z0, xferFields );
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        /* corner at domain logical coord (0, 1, 1) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X1Y0Z0, xferFields)];
        DBGSYNC(xferFields, 1, X1Y0Z0);
        Index_t idx = dx*dy*(dz - 1) + dx*(dy - 1);

        size_t recvCount = xferFields;
        get_yield(src, comBuf, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) += comBuf[fi];
        }
      },
      dash::tasks::in(src),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }
  if (rowNotMax & colNotMax & planeNotMin) {
    int  rank = myRank - domain.tp()*domain.tp() + domain.tp() + 1;
    auto src  = comm.src( rank, X0Y0Z1, xferFields );
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        /* corner at domain logical coord (1, 1, 0) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X0Y0Z1, xferFields)];
        DBGSYNC(xferFields, 1, X0Y0Z1);
        Index_t idx = dx*dy - 1;

        size_t recvCount = xferFields;
        get_yield(src, comBuf, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) += comBuf[fi];
        }
      },
      dash::tasks::in(src),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }
  if (rowNotMax & colNotMax & planeNotMax) {
    int  rank = myRank + domain.tp()*domain.tp() + domain.tp() + 1;
    auto src  = comm.src( rank, X0Y0Z0, xferFields );
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        /* corner at domain logical coord (1, 1, 1) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X0Y0Z0, xferFields)];
        DBGSYNC(xferFields, 1, X0Y0Z0);
        Index_t idx = dx*dy*dz - 1;

        size_t recvCount = xferFields;
        get_yield(src, comBuf, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) += comBuf[fi];
        }
      },
      dash::tasks::in(src),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }
}

void DASHCommSyncPosVel(Domain& domain, DASHComm& comm, int xferFields,
                        Domain_member *fieldData)
{
  if (domain.numRanks() == 1)
    return;

  int myRank;
  bool doRecv = false;
  Index_t maxPlaneComm = xferFields * domain.maxPlaneSize();
  Index_t maxEdgeComm  = xferFields * domain.maxEdgeSize();
  Index_t dx = domain.sizeX() + 1;
  Index_t dy = domain.sizeY() + 1;
  Index_t dz = domain.sizeZ() + 1;
  //MPI_Status status;
  bool rowNotMin, rowNotMax, colNotMin, colNotMax, planeNotMin, planeNotMax;

  /* assume communication to 6 neighbors by default */
  if( domain.rowLoc()   == 0 )               { rowNotMin   = false; }
  if( domain.rowLoc()   == (domain.tp()-1) ) { rowNotMax   = false; }
  if( domain.colLoc()   == 0 )               { colNotMin   = false; }
  if( domain.colLoc()   == (domain.tp()-1) ) { colNotMax   = false; }
  if( domain.planeLoc() == 0 )               { planeNotMin = false; }
  if( domain.planeLoc() == (domain.tp()-1) ) { planeNotMax = false; }

  myRank = dash::myid();
  if (planeNotMin | planeNotMax) {
    /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
    Index_t opCount = dx * dy;

    if (planeNotMin && doRecv) {
      auto src = comm.src(myRank - (domain.tp() * domain.tp()), Z1, xferFields);
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(Z0, xferFields)];

          size_t recvCount = xferFields*opCount;
          get_yield(src, srcAddr, recvCount);

          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<opCount; ++i) {
              (domain.*dest)(i) = srcAddr[i];
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::in(src)
      );
    }
    if (planeNotMax) {
      auto src = comm.src(myRank + (domain.tp() * domain.tp()), Z0, xferFields);
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          // contiguous memory
          srcAddr = &comm.commDataRecv()[comm.offset(Z0, xferFields)];

          size_t recvCount = xferFields*opCount;
          get_yield(src, srcAddr, recvCount);

          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<opCount; ++i) {
              (domain.*dest)(dx*dy*(dz - 1) + i) = srcAddr[i];
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::in(src)
      );
    }
  }
  if (rowNotMin | rowNotMax) {
    /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
    Index_t opCount = dx * dz;

    if (rowNotMin && doRecv) {
      auto src = comm.src(myRank - domain.tp(), Y1, xferFields);
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(Y1, xferFields)];

          size_t recvCount = xferFields*opCount;
          get_yield(src, srcAddr, recvCount);

          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              for (Index_t j=0; j<dx; ++j) {
                (domain.*dest)(i*dx*dy + j) = srcAddr[i*dx + j];
              }
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::in(src)
      );
    }
    if (rowNotMax) {
      auto src = comm.src( myRank + domain.tp(), Y0, xferFields);
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(Y0, xferFields)];

          size_t recvCount = xferFields*opCount;
          get_yield(src, srcAddr, recvCount);

          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              for (Index_t j=0; j<dx; ++j) {
                (domain.*dest)(dx*(dy - 1) + i*dx*dy + j) = srcAddr[i*dx + j];
              }
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::in(src)
      );
    }
  }

  if (colNotMin | colNotMax) {
    /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
    Index_t opCount = dy * dz;

    if (colNotMin && doRecv) {
      auto src = comm.src( myRank - 1, X1, xferFields);
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(X1, xferFields)];

          size_t recvCount = xferFields*opCount;
          get_yield(src, srcAddr, recvCount);

          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              for (Index_t j=0; j<dy; ++j) {
                (domain.*dest)(i*dx*dy + j*dx) = srcAddr[i*dy + j];
              }
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::in(src)
      );
    }
    if (colNotMax) {
      auto src = comm.src( myRank + 1, X0, xferFields );
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(X0, xferFields)];

          size_t recvCount = xferFields*opCount;
          get_yield(src, srcAddr, recvCount);

          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              for (Index_t j=0; j<dy; ++j) {
                (domain.*dest)(dx - 1 + i*dx*dy + j*dx) = srcAddr[i*dy + j];
              }
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::in(src)
      );
    }
  }
  if (rowNotMin && colNotMin && doRecv) {
      int  rank = myRank - domain.tp() - 1;
      auto src  = comm.src( rank, X1Y1, xferFields);
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          srcAddr = &comm.commDataRecv()[comm.offset(X1Y1, xferFields)];

          size_t recvCount = xferFields*dz;
          get_yield(src, srcAddr, recvCount);

          //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              (domain.*dest)(i*dx*dy) = srcAddr[i];
            }
            srcAddr += dz;
          }
        },
        dash::tasks::in(&(domain.*fieldData[0])(0)),
        dash::tasks::in(src)
      );
  }

  if (rowNotMin && planeNotMin && doRecv) {
    int  rank = myRank - domain.tp()*domain.tp() - domain.tp();
    auto src  = comm.src( rank, Y1Z1, xferFields);
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(Y1Z1, xferFields)];

        size_t recvCount = xferFields*dx;
        get_yield(src, srcAddr, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dx; ++i) {
            (domain.*dest)(i) = srcAddr[i];
          }
          srcAddr += dx;
        }
      },
      dash::tasks::in(&(domain.*fieldData[0])(0)),
      dash::tasks::in(src)
    );
  }

  if (colNotMin && planeNotMin && doRecv) {
    int  rank = myRank - domain.tp()*domain.tp() - 1;
    auto src  = comm.src( rank, X1Z1, xferFields);
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X1Z1, xferFields)];

        size_t recvCount = xferFields*dy;
        get_yield(src, srcAddr, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dy; ++i) {
            (domain.*dest)(i*dx) = srcAddr[i];
          }
          srcAddr += dy;
        }
      },
      dash::tasks::in(&(domain.*fieldData[0])(0)),
      dash::tasks::in(src)
    );
  }

  if (rowNotMax && colNotMax) {
    int  rank = myRank + domain.tp() + 1;
    auto src  = comm.src( rank, X0Y0, xferFields);
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X0Y0, xferFields)];

        size_t recvCount = xferFields*dz;
        get_yield(src, srcAddr, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dz; ++i) {
            (domain.*dest)(dx*dy - 1 + i*dx*dy) = srcAddr[i];
          }
          srcAddr += dz;
        }
      },
      dash::tasks::in(&(domain.*fieldData[0])(0)),
      dash::tasks::in(src)
    );
  }

  if (rowNotMax && planeNotMax) {
    int  rank = myRank + domain.tp()*domain.tp() + domain.tp();
    auto src  = comm.src( rank, Y0Z0, xferFields);
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(Y0Z0, xferFields)];

        size_t recvCount = xferFields*dx;
        get_yield(src, srcAddr, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dx; ++i) {
            (domain.*dest)(dx*(dy-1) + dx*dy*(dz-1) + i) = srcAddr[i];
          }
          srcAddr += dx;
        }
      },
      dash::tasks::in(&(domain.*fieldData[0])(0)),
      dash::tasks::in(src)
    );
  }

  if (colNotMax && planeNotMax) {
    int  rank = myRank + domain.tp()*domain.tp() + 1;
    auto src  = comm.src( rank, X0Z0, xferFields);
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X0Z0, xferFields)];

        size_t recvCount = xferFields*dy;
        get_yield(src, srcAddr, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dy; ++i) {
            (domain.*dest)(dx*dy*(dz-1) + dx - 1 + i*dx) = srcAddr[i];
          }
          srcAddr += dy;
        }
      },
      dash::tasks::in(&(domain.*fieldData[0])(0)),
      dash::tasks::in(src)
    );
  }

  if (rowNotMax && colNotMin) {
    int  rank = myRank + domain.tp() - 1;
    auto src  = comm.src( rank, X1Y0, xferFields);
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X1Y0, xferFields)];

        size_t recvCount = xferFields*dz;
        get_yield(src, srcAddr, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dz; ++i) {
            (domain.*dest)(dx*(dy-1) + i*dx*dy) = srcAddr[i];
          }
          srcAddr += dz;
        }
      },
      dash::tasks::in(&(domain.*fieldData[0])(0)),
      dash::tasks::in(src)
    );
  }

  if (rowNotMin && planeNotMax) {
    int  rank = myRank + domain.tp()*domain.tp() - domain.tp();
    auto src  = comm.src( rank, Y1Z0, xferFields);
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(Y1Z0, xferFields)];

        size_t recvCount = xferFields*dx;
        get_yield(src, srcAddr, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dx; ++i) {
            (domain.*dest)(dx*dy*(dz-1) + i) = srcAddr[i];
          }
          srcAddr += dx;
        }
      },
      dash::tasks::in(&(domain.*fieldData[0])(0)),
      dash::tasks::in(src)
    );
  }

  if (colNotMin && planeNotMax) {
    int  rank = myRank + domain.tp()*domain.tp() - 1;
    auto src  = comm.dest( rank, X1Z0, xferFields);
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X1Z0, xferFields)];

        size_t recvCount = xferFields*dy;
        get_yield(src, srcAddr, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dy; ++i) {
            (domain.*dest)(dx*dy*(dz-1) + i*dx) = srcAddr[i];
          }
          srcAddr += dy;
        }
      },
      dash::tasks::in(&(domain.*fieldData[0])(0)),
      dash::tasks::in(src)
    );
  }

  if (rowNotMin && colNotMax && doRecv) {
    int  rank = myRank - domain.tp() + 1;
    auto src  = comm.src( rank, X0Y1, xferFields);
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X0Y1, xferFields)];

        size_t recvCount = xferFields*dz;
        get_yield(src, srcAddr, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dz; ++i) {
            (domain.*dest)(dx - 1 + i*dx*dy) = srcAddr[i];
          }
          srcAddr += dz;
        }
      },
      dash::tasks::in(&(domain.*fieldData[0])(0)),
      dash::tasks::in(src)
    );
  }

  if (rowNotMax && planeNotMin && doRecv) {
    int  rank = myRank - domain.tp()*domain.tp() + domain.tp();
    auto src  = comm.src( rank, Y0Z1, xferFields );
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(Y0Z1, xferFields)];

        size_t recvCount = xferFields*dx;
        get_yield(src, srcAddr, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dx; ++i) {
            (domain.*dest)(dx*(dy - 1) + i) = srcAddr[i];
          }
          srcAddr += dx;
        }
      },
      dash::tasks::in(&(domain.*fieldData[0])(0)),
      dash::tasks::in(src)
    );
  }

  if (colNotMax && planeNotMin && doRecv) {
    int  rank = myRank - domain.tp()*domain.tp() + 1;
    auto src  = comm.dest( rank, X0Z1, xferFields );
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X0Z1, xferFields)];

        size_t recvCount = xferFields*dy;
        get_yield(src, srcAddr, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dy; ++i) {
            (domain.*dest)(dx - 1 + i*dx) = srcAddr[i];
          }
          srcAddr += dy;
        }
      },
      dash::tasks::in(&(domain.*fieldData[0])(0)),
      dash::tasks::in(src)
    );
  }

  if (rowNotMin && colNotMin && planeNotMin && doRecv) {
    int  rank = myRank - domain.tp()*domain.tp() - domain.tp() - 1;
    auto src  = comm.src( rank, X1Y1Z1, xferFields );
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        /* corner at domain logical coord (0, 0, 0) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X1Y1Z1, xferFields)];

        size_t recvCount = xferFields;
        get_yield(src, comBuf, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(0) = comBuf[fi];
        }
      },
      dash::tasks::in(&(domain.*fieldData[0])(0)),
      dash::tasks::in(src)
    );
  }
  if (rowNotMin && colNotMin && planeNotMax) {
    int  rank = myRank + domain.tp()*domain.tp() - domain.tp() - 1;
    auto src  = comm.src( rank, X1Y1Z0, xferFields );
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        /* corner at domain logical coord (0, 0, 1) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X1Y1Z0, xferFields)];
        Index_t idx = dx*dy*(dz - 1);

        size_t recvCount = xferFields;
        get_yield(src, comBuf, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) = comBuf[fi];
        }
      },
      dash::tasks::in(&(domain.*fieldData[0])(0)),
      dash::tasks::in(src)
    );
  }
  if (rowNotMin && colNotMax && planeNotMin && doRecv) {
    int  rank = myRank - domain.tp()*domain.tp() - domain.tp() + 1;
    auto src  = comm.src( rank, X0Y1Z1, xferFields );
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        /* corner at domain logical coord (1, 0, 0) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X0Y1Z1, xferFields)];
        Index_t idx = dx - 1;

        size_t recvCount = xferFields;
        get_yield(src, comBuf, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) = comBuf[fi];
        }
      },
      dash::tasks::in(&(domain.*fieldData[0])(0)),
      dash::tasks::in(src)
    );
  }
  if (rowNotMin && colNotMax && planeNotMax) {
    int  rank = myRank + domain.tp()*domain.tp() - domain.tp() + 1;
    auto src  = comm.src( rank, X0Y1Z0, xferFields );
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        /* corner at domain logical coord (1, 0, 1) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X0Y1Z0, xferFields)];
        Index_t idx = dx*dy*(dz - 1) + (dx - 1);

        size_t recvCount = xferFields;
        get_yield(src, comBuf, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) = comBuf[fi];
        }
      },
      dash::tasks::in(&(domain.*fieldData[0])(0)),
      dash::tasks::in(src)
    );
  }
  if (rowNotMax && colNotMin && planeNotMin && doRecv) {
    int  rank = myRank - domain.tp()*domain.tp() + domain.tp() - 1;
    auto src  = comm.src( rank, X1Y0Z1, xferFields );
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        /* corner at domain logical coord (0, 1, 0) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X1Y0Z1, xferFields)];
        Index_t idx = dx*(dy - 1);

        size_t recvCount = xferFields;
        get_yield(src, comBuf, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) = comBuf[fi];
        }
      },
      dash::tasks::in(&(domain.*fieldData[0])(0)),
      dash::tasks::in(src)
    );
  }
  if (rowNotMax && colNotMin && planeNotMax) {
    int  rank = myRank + domain.tp()*domain.tp() + domain.tp() - 1;
    auto src  = comm.src( rank, X1Y0Z0, xferFields );
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        /* corner at domain logical coord (0, 1, 1) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X1Y0Z0, xferFields)];
        Index_t idx = dx*dy*(dz - 1) + dx*(dy - 1);

        size_t recvCount = xferFields;
        get_yield(src, comBuf, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) = comBuf[fi];
        }
      },
      dash::tasks::in(&(domain.*fieldData[0])(0)),
      dash::tasks::in(src)
    );
  }
  if (rowNotMax && colNotMax && planeNotMin && doRecv) {
    int  rank = myRank - domain.tp()*domain.tp() + domain.tp() + 1;
    auto src  = comm.src( rank, X0Y0Z1, xferFields );
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        /* corner at domain logical coord (1, 1, 0) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X0Y0Z1, xferFields)];
        Index_t idx = dx*dy - 1;

        size_t recvCount = xferFields;
        get_yield(src, comBuf, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) = comBuf[fi];
        }
      },
      dash::tasks::in(&(domain.*fieldData[0])(0)),
      dash::tasks::in(src)
    );
  }
  if (rowNotMax && colNotMax && planeNotMax) {
    int  rank = myRank + domain.tp()*domain.tp() + domain.tp() + 1;
    auto src  = comm.src( rank, X0Y0Z0, xferFields );
    dash::tasks::ASYNC(
      [=, &domain, &comm](){
        /* corner at domain logical coord (1, 1, 1) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X0Y0Z0, xferFields)];
        Index_t idx = dx*dy*dz - 1;

        size_t recvCount = xferFields;
        get_yield(src, comBuf, recvCount);

        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) = comBuf[fi];
        }
      },
      dash::tasks::in(&(domain.*fieldData[0])(0)),
      dash::tasks::in(src)
    );
  }
}

void DASHCommMonoQ(Domain& domain, DASHComm& comm)
{
  if (domain.numRanks() == 1)
    return;

  int myRank;
  Index_t xferFields = 3; /* delv_xi, delv_eta, delv_zeta */
  Domain_member fieldData[3];
  Index_t fieldOffset[3];
  Index_t maxPlaneComm = xferFields * domain.maxPlaneSize();
  Index_t dx = domain.sizeX();
  Index_t dy = domain.sizeY();
  Index_t dz = domain.sizeZ();
  //MPI_Status status;
  bool rowNotMin, rowNotMax, colNotMin, colNotMax, planeNotMin, planeNotMax;

  /* assume communication to 6 neighbors by default */
  if( domain.rowLoc()   == 0 )               { rowNotMin   = false; }
  if( domain.rowLoc()   == (domain.tp()-1) ) { rowNotMax   = false; }
  if( domain.colLoc()   == 0 )               { colNotMin   = false; }
  if( domain.colLoc()   == (domain.tp()-1) ) { colNotMax   = false; }
  if( domain.planeLoc() == 0 )               { planeNotMin = false; }
  if( domain.planeLoc() == (domain.tp()-1) ) { planeNotMax = false; }

  /* point into ghost data area */
  // fieldData[0] = &(domain.delv_xi(domain.numElem()));
  // fieldData[1] = &(domain.delv_eta(domain.numElem()));
  // fieldData[2] = &(domain.delv_zeta(domain.numElem()));
  fieldData[0] = &Domain::delv_xi;
  fieldData[1] = &Domain::delv_eta;
  fieldData[2] = &Domain::delv_zeta;
  fieldOffset[0] = domain.numElem();
  fieldOffset[1] = domain.numElem();
  fieldOffset[2] = domain.numElem();

  myRank = dash::myid();

  if (planeNotMin | planeNotMax) {
    /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
    Index_t opCount = dx * dy;

    if (planeNotMin) {
      auto src = comm.src(myRank - (domain.tp() * domain.tp()), Z1, xferFields);
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(Z1, xferFields)];

          size_t recvCount = xferFields*opCount;
          get_yield(src, srcAddr, recvCount);

          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<opCount; ++i) {
              (domain.*dest)(fieldOffset[fi] + i) = srcAddr[i];
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(src),
        dash::tasks::in(&domain.delv_xi(0))
      );
      for (Index_t fi=0; fi<xferFields; ++fi)
        fieldOffset[fi] += opCount;
    }
    if (planeNotMax) {
      auto src = comm.src(myRank + (domain.tp() * domain.tp()), Z0, xferFields);
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(Z0, xferFields)];

          size_t recvCount = xferFields*opCount;
          auto fut =
            dash::copy_async(
              src,
              src + recvCount,
              srcAddr);
          while(!fut.test()) dash::tasks::yield();

          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<opCount; ++i) {
              (domain.*dest)(fieldOffset[fi] + i) = srcAddr[i];
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(src),
        dash::tasks::in(&domain.delv_xi(0))
      );
      for (Index_t fi=0; fi<xferFields; ++fi)
        fieldOffset[fi] += opCount;
    }
  }

  if (rowNotMin | rowNotMax) {
    /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
    Index_t opCount = dx * dz;

    if (rowNotMin) {
      auto src = comm.src(myRank - domain.tp(), Y1, xferFields);
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(Y1, xferFields)];

          size_t recvCount = xferFields*opCount;
          get_yield(src, srcAddr, recvCount);

          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<opCount; ++i) {
              (domain.*dest)(fieldOffset[fi] + i) = srcAddr[i];
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(src),
        dash::tasks::in(&domain.delv_xi(0))
      );
      for (Index_t fi=0; fi<xferFields; ++fi)
        fieldOffset[fi] += opCount;
    }
    if (rowNotMax) {
      auto src = comm.src(myRank + domain.tp(), Y0, xferFields);
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(Y0, xferFields)];

          size_t recvCount = xferFields*opCount;
          get_yield(src, srcAddr, recvCount);

          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<opCount; ++i) {
              (domain.*dest)(fieldOffset[fi] + i) = srcAddr[i];
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(src),
        dash::tasks::in(&domain.delv_xi(0))
      );
      for (Index_t fi=0; fi<xferFields; ++fi)
        fieldOffset[fi] += opCount;
    }
  }
  if (colNotMin | colNotMax) {
    /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
    Index_t opCount = dy * dz;

    if (colNotMin) {
      auto src = comm.src(myRank - 1, X1, xferFields);
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(X1, xferFields)];

          size_t recvCount = xferFields*opCount;
          get_yield(src, srcAddr, recvCount);

          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<opCount; ++i) {
              (domain.*dest)(fieldOffset[fi] + i) = srcAddr[i];
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(src),
        dash::tasks::in(&domain.delv_xi(0))
      );
      for (Index_t fi=0; fi<xferFields; ++fi)
        fieldOffset[fi] += opCount;
    }
    if (colNotMax) {
      auto src = comm.src(myRank + 1, X0, xferFields);
      dash::tasks::ASYNC(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(X0, xferFields)];

          size_t recvCount = xferFields*opCount;
          get_yield(src, srcAddr, recvCount);

          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<opCount; ++i) {
              (domain.*dest)(fieldOffset[fi] + i) = srcAddr[i];
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(src),
        dash::tasks::in(&domain.delv_xi(0))
      );
    }
  }
}

#endif // DASH_USE_GET
