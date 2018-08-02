
#include <iostream>

#include "lulesh-dash.h"
#include "lulesh-comm-dash.h"
#include "lulesh-comm-dash-onesided.h"

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
      auto dest = comm.dest(myRank - (domain.tp() * domain.tp()), Z1);
      dash::tasks::async(
        [=, &domain, &comm](){
          int sendCount = dx * dy;
          Real_t *destAddr;
          destAddr = &comm.commDataSend()[pmsg * maxPlaneComm];
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<sendCount; ++i) {
              destAddr[i] = (domain.*src)(i);
            }
            destAddr += sendCount;
          }
          destAddr -= xferFields*sendCount;

          auto fut =
            dash::copy_async(
              destAddr,
              destAddr + (xferFields * sendCount),
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
      );
      /*
        MPI_Isend(destAddr, xferFields*sendCount, baseType,
        myRank - domain.tp()*domain.tp(), msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg]);
      */

      ++pmsg;
    }
    if (planeNotMax && doSend) {
      auto dest = comm.dest(myRank + (domain.tp() * domain.tp()), Z0);
      dash::tasks::async(
        [=, &domain, &comm](){
          int sendCount = dx * dy;
          Real_t *destAddr;
          destAddr = &comm.commDataSend()[pmsg * maxPlaneComm];
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<sendCount; ++i) {
              destAddr[i] = (domain.*src)(dx*dy*(dz - 1) + i);
            }
            destAddr += sendCount;
          }
          destAddr -= xferFields*sendCount;

          auto fut =
            dash::copy_async(
              destAddr,
              destAddr + (xferFields * sendCount),
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
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
      auto dest = comm.dest(myRank - domain.tp(), Y1);
      dash::tasks::async(
        [=, &domain, &comm](){
          int sendCount = dx * dz;
          Real_t *destAddr;
          destAddr = &comm.commDataSend()[pmsg * maxPlaneComm];
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              for (Index_t j=0; j<dx; ++j) {
                destAddr[i*dx+j] = (domain.*src)(i*dx*dy + j);
              }
            }
            destAddr += sendCount;
          }
          destAddr -= xferFields*sendCount;

          auto fut =
            dash::copy_async(
              destAddr,
              destAddr + xferFields*sendCount,
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
      );
      /*
        MPI_Isend(destAddr, xferFields*sendCount, baseType,
        myRank - domain.tp(), msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg]);
      */
      ++pmsg;
    }
    if (rowNotMax && doSend) {
      auto dest = comm.dest( myRank + domain.tp(), Y0);
      dash::tasks::async(
        [=, &domain, &comm](){
          int sendCount = dx * dz;
          int myRank = dash::myid();
          Real_t *destAddr;
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
          destAddr -= xferFields*sendCount;

          auto fut =
            dash::copy_async(
              destAddr, destAddr+xferFields*sendCount,
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
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
      auto dest = comm.dest( myRank - 1, X1);
      dash::tasks::async(
        [=, &domain, &comm](){
          int sendCount = dy * dz;
          Real_t *destAddr;
          destAddr = &comm.commDataSend()[pmsg * maxPlaneComm];
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              for (Index_t j=0; j<dy; ++j) {
                destAddr[i*dy + j] = (domain.*src)(i*dx*dy + j*dx);
              }
            }
            destAddr += sendCount;
          }
          destAddr -= xferFields*sendCount;

          auto fut =
            dash::copy_async(
              destAddr, destAddr+xferFields*sendCount,
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
      );

      /*
        MPI_Isend(destAddr, xferFields*sendCount, baseType,
        myRank - 1, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg]);
      */
      ++pmsg;
    }
    if (colNotMax && doSend) {
      auto dest = comm.dest( myRank + 1, X0 );
      dash::tasks::async(
        [=, &domain, &comm](){
          int sendCount = dy * dz;
          Real_t *destAddr;
          destAddr = &comm.commDataSend()[pmsg * maxPlaneComm];
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              for (Index_t j=0; j<dy; ++j) {
                destAddr[i*dy + j] = (domain.*src)(dx - 1 + i*dx*dy + j*dx);
              }
            }
            destAddr += sendCount;
          }
          destAddr -= xferFields*sendCount;

          auto fut =
            dash::copy_async(
              destAddr, destAddr+xferFields*sendCount,
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
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
      int toRank = myRank - domain.tp() - 1;
      auto dest  = comm.dest( toRank, X1Y1);
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *destAddr;
          destAddr = &comm.commDataSend()[pmsg * maxPlaneComm +
                                          emsg * maxEdgeComm];
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              destAddr[i] = (domain.*src)(i*dx*dy);
            }
            destAddr += dz;
          }
          destAddr -= xferFields*dz;

          auto fut =
            dash::copy_async(
              destAddr, destAddr+xferFields*dz,
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
      );

      /*
        MPI_Isend(destAddr, xferFields*dz, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (rowNotMin && planeNotMin) {
      int toRank = myRank - domain.tp()*domain.tp() - domain.tp();
      auto dest = comm.dest( toRank, Y1Z1);
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *destAddr;
          destAddr = &comm.commDataSend()[pmsg * maxPlaneComm +
                                          emsg * maxEdgeComm];
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dx; ++i) {
              destAddr[i] = (domain.*src)(i);
            }
            destAddr += dx;
          }
          destAddr -= xferFields*dx;

          auto fut =
            dash::copy_async(
              destAddr, destAddr+xferFields*dx,
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
      );

      /*
        MPI_Isend(destAddr, xferFields*dx, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (colNotMin && planeNotMin) {
      int toRank = myRank - domain.tp()*domain.tp() - 1;
      auto dest = comm.dest( toRank, X1Z1);
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *destAddr;
          destAddr = &comm.commDataSend()[pmsg * maxPlaneComm +
                                          emsg * maxEdgeComm];
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dy; ++i) {
              destAddr[i] = (domain.*src)(i*dx);
            }
            destAddr += dy;
          }
          destAddr -= xferFields*dy;

          auto fut =
            dash::copy_async(
              destAddr, destAddr+xferFields*dy,
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
      );

      /*
        MPI_Isend(destAddr, xferFields*dy, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (rowNotMax && colNotMax && doSend) {
      int toRank = myRank + domain.tp() + 1;
      auto dest = comm.dest( toRank, X0Y0);
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *destAddr;
          destAddr = &comm.commDataSend()[pmsg * maxPlaneComm +
                                          emsg * maxEdgeComm];
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              destAddr[i] = (domain.*src)(dx*dy - 1 + i*dx*dy);
            }
            destAddr += dz;
          }
          destAddr -= xferFields*dz;

          auto fut =
            dash::copy_async(
              destAddr, destAddr+xferFields*dz,
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
      );
      /*
        MPI_Isend(destAddr, xferFields*dz, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (rowNotMax && planeNotMax && doSend) {
      int toRank = myRank + domain.tp()*domain.tp() + domain.tp();
      auto dest = comm.dest( toRank, Y0Z0);
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *destAddr;
          destAddr = &comm.commDataSend()[pmsg * maxPlaneComm +
                                          emsg * maxEdgeComm];
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dx; ++i) {
              destAddr[i] = (domain.*src)(dx*(dy-1) + dx*dy*(dz-1) + i);
            }
            destAddr += dx;
          }
          destAddr -= xferFields*dx;

          auto fut =
            dash::copy_async(
              destAddr, destAddr+xferFields*dx,
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
      );
      /*
        MPI_Isend(destAddr, xferFields*dx, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (colNotMax && planeNotMax && doSend) {
      int toRank = myRank + domain.tp()*domain.tp() + 1;
      auto dest = comm.dest( toRank, X0Z0);
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *destAddr;
          destAddr = &comm.commDataSend()[pmsg * maxPlaneComm +
                                          emsg * maxEdgeComm];
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dy; ++i) {
              destAddr[i] = (domain.*src)(dx*dy*(dz-1) + dx - 1 + i*dx);
            }
            destAddr += dy;
          }
          destAddr -= xferFields*dy;

          auto fut =
            dash::copy_async(
              destAddr, destAddr+xferFields*dy,
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
      );
      /*
        MPI_Isend(destAddr, xferFields*dy, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (rowNotMax && colNotMin && doSend) {
      int toRank = myRank + domain.tp() - 1;
      auto dest = comm.dest( toRank, X1Y0);
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *destAddr;
          destAddr = &comm.commDataSend()[pmsg * maxPlaneComm +
                                          emsg * maxEdgeComm];
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              destAddr[i] = (domain.*src)(dx*(dy-1) + i*dx*dy);
            }
            destAddr += dz;
          }
          destAddr -= xferFields*dz;

          auto fut =
            dash::copy_async(
              destAddr, destAddr+xferFields*dz,
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
      );
      /*
        MPI_Isend(destAddr, xferFields*dz, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (rowNotMin && planeNotMax && doSend) {
      int toRank = myRank + domain.tp()*domain.tp() - domain.tp();
      auto dest = comm.dest( toRank, Y1Z0);
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *destAddr;
          destAddr = &comm.commDataSend()[pmsg * maxPlaneComm +
                                          emsg * maxEdgeComm];
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dx; ++i) {
              destAddr[i] = (domain.*src)(dx*dy*(dz-1) + i);
            }
            destAddr += dx;
          }
          destAddr -= xferFields*dx;

          auto fut =
            dash::copy_async(
              destAddr, destAddr+xferFields*dx,
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
      );
      /*
        MPI_Isend(destAddr, xferFields*dx, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (colNotMin && planeNotMax && doSend) {
      int toRank = myRank + domain.tp()*domain.tp() - 1;
      auto dest = comm.dest( toRank, X1Z0);
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *destAddr;
          destAddr = &comm.commDataSend()[pmsg * maxPlaneComm +
                                          emsg * maxEdgeComm];
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dy; ++i) {
              destAddr[i] = (domain.*src)(dx*dy*(dz-1) + i*dx);
            }
            destAddr += dy;
          }
          destAddr -= xferFields*dy;

          auto fut =
            dash::copy_async(
              destAddr, destAddr+xferFields*dy,
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
      );
      /*
        MPI_Isend(destAddr, xferFields*dy, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (rowNotMin && colNotMax) {
      int toRank = myRank - domain.tp() + 1;
      auto dest = comm.dest( toRank, X0Y1);
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *destAddr;
          destAddr = &comm.commDataSend()[pmsg * maxPlaneComm +
                                          emsg * maxEdgeComm];
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              destAddr[i] = (domain.*src)(dx - 1 + i*dx*dy);
            }
            destAddr += dz;
          }
          destAddr -= xferFields*dz;

          auto fut =
            dash::copy_async(
              destAddr, destAddr+xferFields*dz,
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
      );
      /*
        MPI_Isend(destAddr, xferFields*dz, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (rowNotMax && planeNotMin) {
      int toRank = myRank - domain.tp()*domain.tp() + domain.tp();
      auto dest = comm.dest( toRank, Y0Z1 );
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *destAddr;
          destAddr = &comm.commDataSend()[pmsg * maxPlaneComm +
                                          emsg * maxEdgeComm];
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dx; ++i) {
              destAddr[i] = (domain.*src)(dx*(dy - 1) + i);
            }
            destAddr += dx;
          }
          destAddr -= xferFields*dx;

          auto fut =
            dash::copy_async(
              destAddr, destAddr+xferFields*dx,
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
      );

      /*
        MPI_Isend(destAddr, xferFields*dx, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (colNotMax && planeNotMin) {
      int toRank = myRank - domain.tp()*domain.tp() + 1;
      auto dest = comm.dest( toRank, X0Z1);
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *destAddr;
          destAddr = &comm.commDataSend()[pmsg * maxPlaneComm +
                                          emsg * maxEdgeComm];
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            for (Index_t i=0; i<dy; ++i) {
              destAddr[i] = (domain.*src)(dx - 1 + i*dx);
            }
            destAddr += dy;
          }
          destAddr -= xferFields*dy;

          auto fut =
            dash::copy_async(
              destAddr, destAddr+xferFields*dy,
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
      );

      /*
        MPI_Isend(destAddr, xferFields*dy, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (rowNotMin && colNotMin && planeNotMin) {
      int toRank = myRank - domain.tp()*domain.tp() - domain.tp() - 1;
      auto dest = comm.dest( toRank, X1Y1Z1);
      dash::tasks::async(
        [=, &domain, &comm](){
          // corner at domain logical coord (0, 0, 0)
          Real_t *comBuf = &comm.commDataSend()[pmsg * maxPlaneComm +
                                                emsg * maxEdgeComm +
                                                cmsg * CACHE_COHERENCE_PAD_REAL];
          for (Index_t fi=0; fi<xferFields; ++fi) {
            comBuf[fi] = (domain.*fieldData[fi])(0);
          }

          auto fut =
            dash::copy_async(
              comBuf, comBuf+xferFields,
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
      );

      /*
        MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]);
      */
      ++cmsg;
    }
    if (rowNotMin && colNotMin && planeNotMax && doSend) {
      int toRank = myRank + domain.tp()*domain.tp() - domain.tp() - 1;
      auto dest = comm.dest( toRank, X1Y1Z0 );
      dash::tasks::async(
        [=, &domain, &comm](){
          // corner at domain logical coord (0, 0, 1)
          Real_t *comBuf = &comm.commDataSend()[pmsg * maxPlaneComm +
                                                emsg * maxEdgeComm +
                                                cmsg * CACHE_COHERENCE_PAD_REAL];
          Index_t idx = dx*dy*(dz - 1);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            comBuf[fi] = (domain.*fieldData[fi])(idx);
          }

          auto fut =
            dash::copy_async(
              comBuf, comBuf+xferFields,
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
      );

      /*
        MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]);
      */
      ++cmsg;
    }
    if (rowNotMin && colNotMax && planeNotMin) {
      int toRank = myRank - domain.tp()*domain.tp() - domain.tp() + 1;
      auto dest = comm.dest( toRank, X0Y1Z1);
      dash::tasks::async(
        [=, &domain, &comm](){
          // corner at domain logical coord (1, 0, 0)
          Real_t *comBuf = &comm.commDataSend()[pmsg * maxPlaneComm +
                                                emsg * maxEdgeComm +
                                                cmsg * CACHE_COHERENCE_PAD_REAL];
          Index_t idx = dx - 1;
          for (Index_t fi=0; fi<xferFields; ++fi) {
            comBuf[fi] = (domain.*fieldData[fi])(idx);
          }

          auto fut =
            dash::copy_async(
              comBuf, comBuf+xferFields,
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
      );

      /*
        MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]);
      */
      ++cmsg;
    }
    if (rowNotMin && colNotMax && planeNotMax && doSend) {
      int toRank = myRank + domain.tp()*domain.tp() - domain.tp() + 1;
      auto dest = comm.dest( toRank, X0Y1Z0 );
      dash::tasks::async(
        [=, &domain, &comm](){
          // corner at domain logical coord (1, 0, 1)
          Real_t *comBuf = &comm.commDataSend()[pmsg * maxPlaneComm +
                                                emsg * maxEdgeComm +
                                                cmsg * CACHE_COHERENCE_PAD_REAL];
          Index_t idx = dx*dy*(dz - 1) + (dx - 1);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            comBuf[fi] = (domain.*fieldData[fi])(idx);
          }

          auto fut =
            dash::copy_async(
              comBuf, comBuf+xferFields,
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
      );

      /*
        MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]);
      */
      ++cmsg;
    }
    if (rowNotMax && colNotMin && planeNotMin) {
      int toRank = myRank - domain.tp()*domain.tp() + domain.tp() - 1;
      auto dest = comm.dest( toRank, X1Y0Z1);
      dash::tasks::async(
        [=, &domain, &comm](){
          // corner at domain logical coord (0, 1, 0)
          Real_t *comBuf = &comm.commDataSend()[pmsg * maxPlaneComm +
                                                emsg * maxEdgeComm +
                                                cmsg * CACHE_COHERENCE_PAD_REAL];
          Index_t idx = dx*(dy - 1);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            comBuf[fi] = (domain.*fieldData[fi])(idx);
          }

          auto fut =
            dash::copy_async(
              comBuf, comBuf+xferFields,
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
      );

      /*
        MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]);
      */
      ++cmsg;
    }
    if (rowNotMax && colNotMin && planeNotMax && doSend) {
      int toRank = myRank + domain.tp()*domain.tp() + domain.tp() - 1;
      auto dest = comm.dest( toRank, X1Y0Z0 );
      dash::tasks::async(
        [=, &domain, &comm](){
          // corner at domain logical coord (0, 1, 1)
          Real_t *comBuf = &comm.commDataSend()[pmsg * maxPlaneComm +
                                                emsg * maxEdgeComm +
                                                cmsg * CACHE_COHERENCE_PAD_REAL];
          Index_t idx = dx*dy*(dz - 1) + dx*(dy - 1);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            comBuf[fi] = (domain.*fieldData[fi])(idx);
          }

          auto fut =
            dash::copy_async(
              comBuf, comBuf+xferFields,
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
      );

      /*
        MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]);
      */
      ++cmsg;
    }
    if (rowNotMax && colNotMax && planeNotMin) {
      int toRank = myRank - domain.tp()*domain.tp() + domain.tp() + 1;
      auto dest = comm.dest( toRank, X0Y0Z1 );
      dash::tasks::async(
        [=, &domain, &comm](){
          // corner at domain logical coord (1, 1, 0)
          Real_t *comBuf = &comm.commDataSend()[pmsg * maxPlaneComm +
                                                emsg * maxEdgeComm +
                                                cmsg * CACHE_COHERENCE_PAD_REAL];
          Index_t idx = dx*dy - 1;
          for (Index_t fi=0; fi<xferFields; ++fi) {
            comBuf[fi] = (domain.*fieldData[fi])(idx);
          }

          auto fut =
            dash::copy_async(
              comBuf, comBuf+xferFields,
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
      );

      /*
        MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]);
      */
      ++cmsg;
    }
    if (rowNotMax && colNotMax && planeNotMax && doSend) {
      int toRank = myRank + domain.tp()*domain.tp() + domain.tp() + 1;
      auto dest = comm.dest( toRank, X0Y0Z0 );
      dash::tasks::async(
        [=, &domain, &comm](){
          // corner at domain logical coord (1, 1, 1)
          Real_t *comBuf = &comm.commDataSend()[pmsg * maxPlaneComm +
                                                emsg * maxEdgeComm +
                                                cmsg * CACHE_COHERENCE_PAD_REAL];
          Index_t idx = dx*dy*dz - 1;
          for (Index_t fi=0; fi<xferFields; ++fi) {
            comBuf[fi] = (domain.*fieldData[fi])(idx);
          }

          auto fut =
            dash::copy_async(
              comBuf, comBuf+xferFields,
              dest);
          while(!fut.test()) dash::tasks::yield();
        },
        [=, &domain](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member src = fieldData[fi];
            *deps = dash::tasks::in(&(domain.*src)(0));
          }
          *deps = dash::tasks::out(dest);
        }
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
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(Z0)];
          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          DBGSYNC(xferFields, opCount, Z0);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<opCount; ++i) {
              (domain.*dest)(i) += srcAddr[i];
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(&comm.commDataRecv()[comm.offset(Z0)]),
        // TODO: dummy IN dependencies, replace them with CONCURRENT!!
        dash::tasks::in(&(domain.*fieldData[0])(0))
      );
    }
    if (planeNotMax) {
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(Z1)];
          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          DBGSYNC(xferFields, opCount, Z1);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<opCount; ++i) {
              (domain.*dest)(dx*dy*(dz - 1) + i) += srcAddr[i];
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(&comm.commDataRecv()[comm.offset(Z1)]),
        // TODO: dummy IN dependencies, replace them with CONCURRENT!!
        dash::tasks::in(&(domain.*fieldData[0])(0))
      );
    }
  }

  if (rowNotMin | rowNotMax) {
    /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
    Index_t opCount = dx * dz;

    if (rowNotMin) {
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(Y0)];
          DBGSYNC(xferFields, opCount, Y0);
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
        dash::tasks::in(&comm.commDataRecv()[comm.offset(Y0)]),
        // TODO: dummy IN dependencies, replace them with CONCURRENT!!
        dash::tasks::in(&(domain.*fieldData[0])(0))
      );
    }
    if (rowNotMax) {
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(Y1)];
          DBGSYNC(xferFields, opCount, Y1);
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
        dash::tasks::in(&comm.commDataRecv()[comm.offset(Y1)]),
        // TODO: dummy IN dependencies, replace them with CONCURRENT!!
        dash::tasks::in(&(domain.*fieldData[0])(0))
      );
    }
  }
  if (colNotMin | colNotMax) {
    /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
    Index_t opCount = dy * dz;

    if (colNotMin) {
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(X0)];
          DBGSYNC(xferFields, opCount, X0);
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
        dash::tasks::in(&comm.commDataRecv()[comm.offset(X0)]),
        // TODO: dummy IN dependencies, replace them with CONCURRENT!!
        dash::tasks::in(&(domain.*fieldData[0])(0))
      );
    }
    if (colNotMax) {
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(X1)];
          DBGSYNC(xferFields, opCount, X1);
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
        dash::tasks::in(&comm.commDataRecv()[comm.offset(X1)]),
        // TODO: dummy IN dependencies, replace them with CONCURRENT!!
        dash::tasks::in(&(domain.*fieldData[0])(0))
      );
    }
  }

  if (rowNotMin & colNotMin) {
    dash::tasks::async(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X0Y0)];
        DBGSYNC(xferFields, dz, X0Y0);
        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dz; ++i) {
            (domain.*dest)(i*dx*dy) += srcAddr[i];
          }
          srcAddr += dz;
        }
      },
      dash::tasks::in(&comm.commDataRecv()[comm.offset(X0Y0)]),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }

  if (rowNotMin & planeNotMin) {
    dash::tasks::async(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(Y0Z0)];
        DBGSYNC(xferFields, dx, Y0Z0);
        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dx; ++i) {
            (domain.*dest)(i) += srcAddr[i];
          }
          srcAddr += dx;
        }
      },
      dash::tasks::in(&comm.commDataRecv()[comm.offset(Y0Z0)]),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }

  if (colNotMin & planeNotMin) {
    dash::tasks::async(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X0Z0)];
        DBGSYNC(xferFields, dy, X0Z0);
        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dy; ++i) {
            (domain.*dest)(i*dx) += srcAddr[i];
          }
          srcAddr += dy;
        }
      },
      dash::tasks::in(&comm.commDataRecv()[comm.offset(X0Z0)]),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }

  if (rowNotMax & colNotMax) {
    dash::tasks::async(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X1Y1)];
        DBGSYNC(xferFields, dz, X1Y1);
        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dz; ++i) {
            (domain.*dest)(dx*dy - 1 + i*dx*dy) += srcAddr[i];
          }
          srcAddr += dz;
        }
      },
      dash::tasks::in(&comm.commDataRecv()[comm.offset(X1Y1)]),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }

  if (rowNotMax & planeNotMax) {
    dash::tasks::async(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(Y1Z1)];
        DBGSYNC(xferFields, dx, Y1Z1);
        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dx; ++i) {
            (domain.*dest)(dx*(dy-1) + dx*dy*(dz-1) + i) += srcAddr[i];
          }
          srcAddr += dx;
        }
      },
      dash::tasks::in(&comm.commDataRecv()[comm.offset(Y1Z1)]),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }

  if (colNotMax & planeNotMax) {
    dash::tasks::async(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X1Z1)];
        DBGSYNC(xferFields, dy, X1Z1);
        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dy; ++i) {
            (domain.*dest)(dx*dy*(dz-1) + dx - 1 + i*dx) += srcAddr[i];
          }
          srcAddr += dy;
        }
      },
      dash::tasks::in(&comm.commDataRecv()[comm.offset(X1Z1)]),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }

  if (rowNotMax & colNotMin) {
    dash::tasks::async(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X0Y1)];
        DBGSYNC(xferFields, dz, X0Y1);
        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dz; ++i) {
            (domain.*dest)(dx*(dy-1) + i*dx*dy) += srcAddr[i];
          }
          srcAddr += dz;
        }
      },
      dash::tasks::in(&comm.commDataRecv()[comm.offset(X0Y1)]),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }

  if (rowNotMin & planeNotMax) {
    dash::tasks::async(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(Y0Z1)];
        DBGSYNC(xferFields, dx, Y0Z1);
        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dx; ++i) {
            (domain.*dest)(dx*dy*(dz-1) + i) += srcAddr[i];
          }
          srcAddr += dx;
        }
      },
      dash::tasks::in(&comm.commDataRecv()[comm.offset(Y0Z1)]),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }

  if (colNotMin & planeNotMax) {
    dash::tasks::async(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X0Z1)];
        DBGSYNC(xferFields, dy, X0Z1);
        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dy; ++i) {
            (domain.*dest)(dx*dy*(dz-1) + i*dx) += srcAddr[i];
          }
          srcAddr += dy;
        }
      },
      dash::tasks::in(&comm.commDataRecv()[comm.offset(X0Z1)]),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }

  if (rowNotMin & colNotMax) {
    dash::tasks::async(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X1Y0)];
        DBGSYNC(xferFields, dz, X1Y0);
        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dz; ++i) {
            (domain.*dest)(dx - 1 + i*dx*dy) += srcAddr[i];
          }
          srcAddr += dz;
        }
      },
      dash::tasks::in(&comm.commDataRecv()[comm.offset(X1Y0)]),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }

  if (rowNotMax & planeNotMin) {
    dash::tasks::async(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(Y1Z0)];
        DBGSYNC(xferFields, dx, Y1Z0);
        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dx; ++i) {
            (domain.*dest)(dx*(dy - 1) + i) += srcAddr[i];
          }
          srcAddr += dx;
        }
      },
      dash::tasks::in(&comm.commDataRecv()[comm.offset(Y1Z0)]),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }

  if (colNotMax & planeNotMin) {
    dash::tasks::async(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X1Z0)];
        DBGSYNC(xferFields, dy, X1Z0);
        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dy; ++i) {
            (domain.*dest)(dx - 1 + i*dx) += srcAddr[i];
          }
          srcAddr += dy;
        }
      },
      dash::tasks::in(&comm.commDataRecv()[comm.offset(X1Z0)]),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }

  if (rowNotMin & colNotMin & planeNotMin) {
    dash::tasks::async(
      [=, &domain, &comm](){
        /* corner at domain logical coord (0, 0, 0) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X0Y0Z0)];
        DBGSYNC(xferFields, 1, X0Y0Z0);
        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(0) += comBuf[fi];
        }
      },
      dash::tasks::in(&comm.commDataRecv()[comm.offset(X0Y0Z0)]),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }
  if (rowNotMin & colNotMin & planeNotMax) {
    dash::tasks::async(
      [=, &domain, &comm](){
        /* corner at domain logical coord (0, 0, 1) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X0Y0Z1)];
        DBGSYNC(xferFields, 1, X0Y0Z1);
        Index_t idx = dx*dy*(dz - 1);
        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) += comBuf[fi];
        }
      },
      dash::tasks::in(&comm.commDataRecv()[comm.offset(X0Y0Z1)]),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }
  if (rowNotMin & colNotMax & planeNotMin) {
    dash::tasks::async(
      [=, &domain, &comm](){
        /* corner at domain logical coord (1, 0, 0) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X1Y0Z0)];
        DBGSYNC(xferFields, 1, X1Y0Z0);
        Index_t idx = dx - 1;
        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) += comBuf[fi];
        }
      },
      dash::tasks::in(&comm.commDataRecv()[comm.offset(X1Y0Z0)]),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }
  if (rowNotMin & colNotMax & planeNotMax) {
    dash::tasks::async(
      [=, &domain, &comm](){
        /* corner at domain logical coord (1, 0, 1) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X1Y0Z1)];
        DBGSYNC(xferFields, 1, X1Y0Z1);
        Index_t idx = dx*dy*(dz - 1) + (dx - 1);
        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) += comBuf[fi];
        }
      },
      dash::tasks::in(&comm.commDataRecv()[comm.offset(X1Y0Z1)]),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }
  if (rowNotMax & colNotMin & planeNotMin) {
    dash::tasks::async(
      [=, &domain, &comm](){
        /* corner at domain logical coord (0, 1, 0) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X0Y1Z0)];
        DBGSYNC(xferFields, 1, X0Y1Z0);
        Index_t idx = dx*(dy - 1);
        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) += comBuf[fi];
        }
      },
      dash::tasks::in(&comm.commDataRecv()[comm.offset(X0Y1Z0)]),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }
  if (rowNotMax & colNotMin & planeNotMax) {
    dash::tasks::async(
      [=, &domain, &comm](){
        /* corner at domain logical coord (0, 1, 1) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X0Y1Z1)];
        DBGSYNC(xferFields, 1, X0Y1Z1);
        Index_t idx = dx*dy*(dz - 1) + dx*(dy - 1);
        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) += comBuf[fi];
        }
      },
      dash::tasks::in(&comm.commDataRecv()[comm.offset(X0Y1Z1)]),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }
  if (rowNotMax & colNotMax & planeNotMin) {
    dash::tasks::async(
      [=, &domain, &comm](){
        /* corner at domain logical coord (1, 1, 0) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X1Y1Z0)];
        DBGSYNC(xferFields, 1, X1Y1Z0);
        Index_t idx = dx*dy - 1;
        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) += comBuf[fi];
        }
      },
      dash::tasks::in(&comm.commDataRecv()[comm.offset(X1Y1Z0)]),
      // TODO: dummy IN dependencies, replace them with CONCURRENT!!
      dash::tasks::in(&(domain.*fieldData[0])(0))
    );
  }
  if (rowNotMax & colNotMax & planeNotMax) {
    dash::tasks::async(
      [=, &domain, &comm](){
        /* corner at domain logical coord (1, 1, 1) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X1Y1Z1)];
        DBGSYNC(xferFields, 1, X1Y1Z1);
        Index_t idx = dx*dy*dz - 1;
        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) += comBuf[fi];
        }
      },
      dash::tasks::in(&comm.commDataRecv()[comm.offset(X1Y1Z1)]),
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
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(Z0)];
          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<opCount; ++i) {
              (domain.*dest)(i) = srcAddr[i];
            }
            srcAddr += opCount;
          }
        },
        [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
            *deps = dash::tasks::in(&(domain.*dest)(0));
          }
          *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(Z0)]);
        }
      );
    }
    if (planeNotMax) {
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          // contiguous memory
          srcAddr = &comm.commDataRecv()[comm.offset(Z1)];
          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<opCount; ++i) {
              (domain.*dest)(dx*dy*(dz - 1) + i) = srcAddr[i];
            }
            srcAddr += opCount;
          }
        },
        [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
            *deps = dash::tasks::in(&(domain.*dest)(0));
          }
          *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(Z1)]);
        }
      );
    }
  }
  if (rowNotMin | rowNotMax) {
    /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
    Index_t opCount = dx * dz;

    if (rowNotMin && doRecv) {
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(Y0)];
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
        [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
            *deps = dash::tasks::in(&(domain.*dest)(0));
          }
          *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(Y0)]);
        }
      );
    }
    if (rowNotMax) {
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(Y1)];
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
        [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
            *deps = dash::tasks::in(&(domain.*dest)(0));
          }
          *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(Y1)]);
        }
      );
    }
  }

  if (colNotMin | colNotMax) {
    /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
    Index_t opCount = dy * dz;

    if (colNotMin && doRecv) {
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(X0)];
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
        [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
            *deps = dash::tasks::in(&(domain.*dest)(0));
          }
          *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(X0)]);
        }
      );
    }
    if (colNotMax) {
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(X1)];
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
        [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
            *deps = dash::tasks::in(&(domain.*dest)(0));
          }
          *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(X1)]);
        }
      );
    }
  }
  if (rowNotMin && colNotMin && doRecv) {
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          srcAddr = &comm.commDataRecv()[comm.offset(X0Y0)];
          //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<dz; ++i) {
              (domain.*dest)(i*dx*dy) = srcAddr[i];
            }
            srcAddr += dz;
          }
        },
        [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
            *deps = dash::tasks::in(&(domain.*dest)(0));
          }
          *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(X0Y0)]);
        }
      );
  }

  if (rowNotMin && planeNotMin && doRecv) {
    dash::tasks::async(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X0Y0)];
        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dx; ++i) {
            (domain.*dest)(i) = srcAddr[i];
          }
          srcAddr += dx;
        }
      },
      [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
          *deps = dash::tasks::in(&(domain.*dest)(0));
        }
        *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(X0Y0)]);
      }
    );
  }

  if (colNotMin && planeNotMin && doRecv) {
    dash::tasks::async(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X0Z0)];
        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dy; ++i) {
            (domain.*dest)(i*dx) = srcAddr[i];
          }
          srcAddr += dy;
        }
      },
      [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
          *deps = dash::tasks::in(&(domain.*dest)(0));
        }
        *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(X0Z0)]);
      }
    );
  }

  if (rowNotMax && colNotMax) {
    dash::tasks::async(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X1Y1)];
        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dz; ++i) {
            (domain.*dest)(dx*dy - 1 + i*dx*dy) = srcAddr[i];
          }
          srcAddr += dz;
        }
      },
      [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
          *deps = dash::tasks::in(&(domain.*dest)(0));
        }
        *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(X1Y1)]);
      }
    );
  }

  if (rowNotMax && planeNotMax) {
    dash::tasks::async(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(Y1Z1)];
        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dx; ++i) {
            (domain.*dest)(dx*(dy-1) + dx*dy*(dz-1) + i) = srcAddr[i];
          }
          srcAddr += dx;
        }
      },
      [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
          *deps = dash::tasks::in(&(domain.*dest)(0));
        }
        *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(Y1Z1)]);
      }
    );
  }

  if (colNotMax && planeNotMax) {
    dash::tasks::async(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X1Z1)];
        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dy; ++i) {
            (domain.*dest)(dx*dy*(dz-1) + dx - 1 + i*dx) = srcAddr[i];
          }
          srcAddr += dy;
        }
      },
      [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
          *deps = dash::tasks::in(&(domain.*dest)(0));
        }
        *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(X1Z1)]);
      }
    );
  }

  if (rowNotMax && colNotMin) {
    dash::tasks::async(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X0Y1)];
        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dz; ++i) {
            (domain.*dest)(dx*(dy-1) + i*dx*dy) = srcAddr[i];
          }
          srcAddr += dz;
        }
      },
      [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
          *deps = dash::tasks::in(&(domain.*dest)(0));
        }
        *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(X0Y1)]);
      }
    );
  }

  if (rowNotMin && planeNotMax) {
    dash::tasks::async(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(Y0Z1)];
        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dx; ++i) {
            (domain.*dest)(dx*dy*(dz-1) + i) = srcAddr[i];
          }
          srcAddr += dx;
        }
      },
      [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
          *deps = dash::tasks::in(&(domain.*dest)(0));
        }
        *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(Y0Z1)]);
      }
    );
  }

  if (colNotMin && planeNotMax) {
    dash::tasks::async(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X0Z1)];
        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dy; ++i) {
            (domain.*dest)(dx*dy*(dz-1) + i*dx) = srcAddr[i];
          }
          srcAddr += dy;
        }
      },
      [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
          *deps = dash::tasks::in(&(domain.*dest)(0));
        }
        *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(X0Z1)]);
      }
    );
  }

  if (rowNotMin && colNotMax && doRecv) {
    dash::tasks::async(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X1Y0)];
        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dz; ++i) {
            (domain.*dest)(dx - 1 + i*dx*dy) = srcAddr[i];
          }
          srcAddr += dz;
        }
      },
      [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
          *deps = dash::tasks::in(&(domain.*dest)(0));
        }
        *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(X1Y0)]);
      }
    );
  }

  if (rowNotMax && planeNotMin && doRecv) {
    dash::tasks::async(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(Y1Z0)];
        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dx; ++i) {
            (domain.*dest)(dx*(dy - 1) + i) = srcAddr[i];
          }
          srcAddr += dx;
        }
      },
      [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
          *deps = dash::tasks::in(&(domain.*dest)(0));
        }
        *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(Y1Z0)]);
      }
    );
  }

  if (colNotMax && planeNotMin && doRecv) {
    dash::tasks::async(
      [=, &domain, &comm](){
        Real_t *srcAddr;
        srcAddr = &comm.commDataRecv()[comm.offset(X1Z0)];
        //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          for (Index_t i=0; i<dy; ++i) {
            (domain.*dest)(dx - 1 + i*dx) = srcAddr[i];
          }
          srcAddr += dy;
        }
      },
      [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
          *deps = dash::tasks::in(&(domain.*dest)(0));
        }
        *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(X1Z0)]);
      }
    );
  }

  if (rowNotMin && colNotMin && planeNotMin && doRecv) {
    dash::tasks::async(
      [=, &domain, &comm](){
        /* corner at domain logical coord (0, 0, 0) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X0Y0Z0)];
        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(0) = comBuf[fi];
        }
      },
      [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
          *deps = dash::tasks::in(&(domain.*dest)(0));
        }
        *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(X0Y0Z0)]);
      }
    );
  }
  if (rowNotMin && colNotMin && planeNotMax) {
    dash::tasks::async(
      [=, &domain, &comm](){
        /* corner at domain logical coord (0, 0, 1) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X0Y0Z1)];
        Index_t idx = dx*dy*(dz - 1);
        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) = comBuf[fi];
        }
      },
      [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
          *deps = dash::tasks::in(&(domain.*dest)(0));
        }
        *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(X0Y0Z1)]);
      }
    );
  }
  if (rowNotMin && colNotMax && planeNotMin && doRecv) {
    dash::tasks::async(
      [=, &domain, &comm](){
        /* corner at domain logical coord (1, 0, 0) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X1Y0Z0)];
        Index_t idx = dx - 1;
        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) = comBuf[fi];
        }
      },
      [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
          *deps = dash::tasks::in(&(domain.*dest)(0));
        }
        *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(X1Y0Z0)]);
      }
    );
  }
  if (rowNotMin && colNotMax && planeNotMax) {
    dash::tasks::async(
      [=, &domain, &comm](){
        /* corner at domain logical coord (1, 0, 1) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X1Y0Z1)];
        Index_t idx = dx*dy*(dz - 1) + (dx - 1);
        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) = comBuf[fi];
        }
      },
      [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
          *deps = dash::tasks::in(&(domain.*dest)(0));
        }
        *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(X1Y0Z1)]);
      }
    );
  }
  if (rowNotMax && colNotMin && planeNotMin && doRecv) {
    dash::tasks::async(
      [=, &domain, &comm](){
        /* corner at domain logical coord (0, 1, 0) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X0Y1Z0)];
        Index_t idx = dx*(dy - 1);
        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) = comBuf[fi];
        }
      },
      [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
          *deps = dash::tasks::in(&(domain.*dest)(0));
        }
        *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(X0Y1Z0)]);
      }
    );
  }
  if (rowNotMax && colNotMin && planeNotMax) {
    dash::tasks::async(
      [=, &domain, &comm](){
        /* corner at domain logical coord (0, 1, 1) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X0Y1Z1)];
        Index_t idx = dx*dy*(dz - 1) + dx*(dy - 1);
        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) = comBuf[fi];
        }
      },
      [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
          *deps = dash::tasks::in(&(domain.*dest)(0));
        }
        *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(X0Y1Z1)]);
      }
    );
  }
  if (rowNotMax && colNotMax && planeNotMin && doRecv) {
    dash::tasks::async(
      [=, &domain, &comm](){
        /* corner at domain logical coord (1, 1, 0) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X1Y1Z0)];
        Index_t idx = dx*dy - 1;
        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) = comBuf[fi];
        }
      },
      [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
          *deps = dash::tasks::in(&(domain.*dest)(0));
        }
        *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(X1Y1Z0)]);
      }
    );
  }
  if (rowNotMax && colNotMax && planeNotMax) {
    dash::tasks::async(
      [=, &domain, &comm](){
        /* corner at domain logical coord (1, 1, 1) */
        Real_t *comBuf = &comm.commDataRecv()[comm.offset(X1Y1Z1)];
        Index_t idx = dx*dy*dz - 1;
        //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
        for (Index_t fi=0; fi<xferFields; ++fi) {
          (domain.*fieldData[fi])(idx) = comBuf[fi];
        }
      },
      [=, &domain, &comm](dash::tasks::DependencyVectorInserter deps){
        for (Index_t fi=0; fi<xferFields; ++fi) {
          Domain_member dest = fieldData[fi];
          // TODO: these are dummy IN dependencies, replace them with CONCURRENT!!
          *deps = dash::tasks::in(&(domain.*dest)(0));
        }
        *deps = dash::tasks::in(&comm.commDataRecv()[comm.offset(X1Y1Z1)]);
      }
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
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(Z0)];
          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<opCount; ++i) {
              (domain.*dest)(fieldOffset[fi] + i) = srcAddr[i];
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(&comm.commDataRecv()[comm.offset(Z0)]),
        dash::tasks::in(&domain.delv_xi(0))
      );
      for (Index_t fi=0; fi<xferFields; ++fi)
        fieldOffset[fi] += opCount;
    }
    if (planeNotMax) {
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(Z1)];
          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<opCount; ++i) {
              (domain.*dest)(fieldOffset[fi] + i) = srcAddr[i];
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(&comm.commDataRecv()[comm.offset(Z0)]),
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
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(Y0)];
          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<opCount; ++i) {
              (domain.*dest)(fieldOffset[fi] + i) = srcAddr[i];
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(&comm.commDataRecv()[comm.offset(Z0)]),
        dash::tasks::in(&domain.delv_xi(0))
      );
      for (Index_t fi=0; fi<xferFields; ++fi)
        fieldOffset[fi] += opCount;
    }
    if (rowNotMax) {
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(Y1)];
          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<opCount; ++i) {
              (domain.*dest)(fieldOffset[fi] + i) = srcAddr[i];
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(&comm.commDataRecv()[comm.offset(Y1)]),
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
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(X0)];
          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<opCount; ++i) {
              (domain.*dest)(fieldOffset[fi] + i) = srcAddr[i];
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(&comm.commDataRecv()[comm.offset(X0)]),
        dash::tasks::in(&domain.delv_xi(0))
      );
      for (Index_t fi=0; fi<xferFields; ++fi)
        fieldOffset[fi] += opCount;
    }
    if (colNotMax) {
      dash::tasks::async(
        [=, &domain, &comm](){
          Real_t *srcAddr;
          /* contiguous memory */
          srcAddr = &comm.commDataRecv()[comm.offset(X1)];
          //MPI_Wait(&comm.recvRequest[pmsg], &status);
          for (Index_t fi=0; fi<xferFields; ++fi) {
            Domain_member dest = fieldData[fi];
            for (Index_t i=0; i<opCount; ++i) {
              (domain.*dest)(fieldOffset[fi] + i) = srcAddr[i];
            }
            srcAddr += opCount;
          }
        },
        dash::tasks::in(&comm.commDataRecv()[comm.offset(X1)]),
        dash::tasks::in(&domain.delv_xi(0))
      );
    }
  }
}
