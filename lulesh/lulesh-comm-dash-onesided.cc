
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
  // assume communication to 6 neighbors by default
  rowNotMin = rowNotMax = colNotMin = colNotMax = planeNotMin = planeNotMax = true;

  Index_t maxPlaneComm = xferFields * domain.maxPlaneSize();
  Index_t maxEdgeComm = xferFields * domain.maxEdgeSize();

  Index_t pmsg = 0; // plane comm msg
  Index_t emsg = 0; // edge comm msg
  Index_t cmsg = 0; // corner comm msg


  if( domain.rowLoc()   == 0 )               { rowNotMin   = false; }
  if( domain.rowLoc()   == (domain.tp()-1) ) { rowNotMax   = false; }
  if( domain.colLoc()   == 0 )               { colNotMin   = false; }
  if( domain.colLoc()   == (domain.tp()-1) ) { colNotMax   = false; }
  if( domain.planeLoc() == 0 )               { planeNotMin = false; }
  if( domain.planeLoc() == (domain.tp()-1) ) { planeNotMax = false; }

  int myRank = dash::myid();

  if( planeNotMin | planeNotMax ) {
    // ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE
    int sendCount = dx * dy;

    if (planeNotMin) {
      destAddr = &comm.commDataSend()[pmsg * maxPlaneComm];
      for (Index_t fi=0; fi<xferFields; ++fi) {
        Domain_member src = fieldData[fi];
        for (Index_t i=0; i<sendCount; ++i) {
          destAddr[i] = (domain.*src)(i);
        }
        destAddr += sendCount;
      }
      destAddr -= xferFields*sendCount;

      comm.sendRequest[pmsg] =
        dash::copy_async(
          destAddr,
          destAddr + (xferFields * sendCount),
          comm.dest(myRank - (domain.tp() * domain.tp()), Z1));

      /*
        MPI_Isend(destAddr, xferFields*sendCount, baseType,
        myRank - domain.tp()*domain.tp(), msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg]);
      */

      ++pmsg;
    }
    if (planeNotMax && doSend) {
      destAddr = &comm.commDataSend()[pmsg * maxPlaneComm];
      for (Index_t fi=0; fi<xferFields; ++fi) {
        Domain_member src = fieldData[fi];
        for (Index_t i=0; i<sendCount; ++i) {
          destAddr[i] = (domain.*src)(dx*dy*(dz - 1) + i);
        }
        destAddr += sendCount;
      }
      destAddr -= xferFields*sendCount;

      comm.sendRequest[pmsg] =
        dash::copy_async(
          destAddr,
          destAddr + (xferFields * sendCount),
          comm.dest(myRank + (domain.tp() * domain.tp()), Z0));
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
    int sendCount = dx * dz;

    if (rowNotMin) {
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

      comm.sendRequest[pmsg] =
        dash::copy_async(
          destAddr,
          destAddr + xferFields*sendCount,
          comm.dest(myRank - domain.tp(), Y1));
      /*
        MPI_Isend(destAddr, xferFields*sendCount, baseType,
        myRank - domain.tp(), msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg]);
      */
      ++pmsg;
    }
    if (rowNotMax && doSend) {
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

      comm.sendRequest[pmsg] =
        dash::copy_async(
          destAddr, destAddr+xferFields*sendCount,
          comm.dest( myRank + domain.tp(), Y0 ));

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
    int sendCount = dy * dz;

    if (colNotMin) {
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

      comm.sendRequest[pmsg] =
        dash::copy_async(
          destAddr, destAddr+xferFields*sendCount,
          comm.dest( myRank - 1, X1));
      /*
        MPI_Isend(destAddr, xferFields*sendCount, baseType,
        myRank - 1, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg]);
      */
      ++pmsg;
    }
    if (colNotMax && doSend) {
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

      comm.sendRequest[pmsg] =
        dash::copy_async(
          destAddr, destAddr+xferFields*sendCount,
          comm.dest( myRank + 1, X0 ));
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

      comm.sendRequest[pmsg + emsg] =
        dash::copy_async(
          destAddr, destAddr+xferFields*dz,
          comm.dest( toRank, X1Y1) );

      /*
        MPI_Isend(destAddr, xferFields*dz, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (rowNotMin && planeNotMin) {
      int toRank = myRank - domain.tp()*domain.tp() - domain.tp();
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

      comm.sendRequest[pmsg + emsg] =
        dash::copy_async(
          destAddr, destAddr+xferFields*dx,
          comm.dest( toRank, Y1Z1) );
      /*
        MPI_Isend(destAddr, xferFields*dx, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (colNotMin && planeNotMin) {
      int toRank = myRank - domain.tp()*domain.tp() - 1;
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

      comm.sendRequest[pmsg + emsg] =
        dash::copy_async(
          destAddr, destAddr+xferFields*dy,
          comm.dest( toRank, X1Z1) );
      /*
        MPI_Isend(destAddr, xferFields*dy, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (rowNotMax && colNotMax && doSend) {
      int toRank = myRank + domain.tp() + 1;
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

      comm.sendRequest[pmsg + emsg] =
        dash::copy_async(
          destAddr, destAddr+xferFields*dz,
          comm.dest( toRank, X0Y0) );
      /*
        MPI_Isend(destAddr, xferFields*dz, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (rowNotMax && planeNotMax && doSend) {
      int toRank = myRank + domain.tp()*domain.tp() + domain.tp();
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

      comm.sendRequest[pmsg + emsg] =
        dash::copy_async(
          destAddr, destAddr+xferFields*dx,
          comm.dest( toRank, Y0Z0) );
      /*
        MPI_Isend(destAddr, xferFields*dx, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (colNotMax && planeNotMax && doSend) {
      int toRank = myRank + domain.tp()*domain.tp() + 1;
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

      comm.sendRequest[pmsg + emsg] =
        dash::copy_async(
          destAddr, destAddr+xferFields*dy,
          comm.dest( toRank, X0Z0) );
      /*
        MPI_Isend(destAddr, xferFields*dy, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (rowNotMax && colNotMin && doSend) {
      int toRank = myRank + domain.tp() - 1;
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

      comm.sendRequest[pmsg + emsg] =
        dash::copy_async(
          destAddr, destAddr+xferFields*dz,
          comm.dest( toRank, X1Y0) );
      /*
        MPI_Isend(destAddr, xferFields*dz, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (rowNotMin && planeNotMax && doSend) {
      int toRank = myRank + domain.tp()*domain.tp() - domain.tp();
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

      comm.sendRequest[pmsg + emsg] =
        dash::copy_async(
          destAddr, destAddr+xferFields*dx,
          comm.dest( toRank, Y1Z0) );
      /*
        MPI_Isend(destAddr, xferFields*dx, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (colNotMin && planeNotMax && doSend) {
      int toRank = myRank + domain.tp()*domain.tp() - 1;
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

      comm.sendRequest[pmsg + emsg] =
        dash::copy_async(
          destAddr, destAddr+xferFields*dy,
          comm.dest( toRank, X1Z0) );
      /*
        MPI_Isend(destAddr, xferFields*dy, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (rowNotMin && colNotMax) {
      int toRank = myRank - domain.tp() + 1;
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

      comm.sendRequest[pmsg + emsg] =
        dash::copy_async(
          destAddr, destAddr+xferFields*dz,
          comm.dest( toRank, X0Y1) );
      /*
        MPI_Isend(destAddr, xferFields*dz, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (rowNotMax && planeNotMin) {
      int toRank = myRank - domain.tp()*domain.tp() + domain.tp();
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

      comm.sendRequest[pmsg + emsg] =
        dash::copy_async(
          destAddr, destAddr+xferFields*dx,
          comm.dest( toRank, Y0Z1 ) );

      /*
        MPI_Isend(destAddr, xferFields*dx, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (colNotMax && planeNotMin) {
      int toRank = myRank - domain.tp()*domain.tp() + 1;
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

      comm.sendRequest[pmsg + emsg] =
        dash::copy_async(
          destAddr, destAddr+xferFields*dy,
          comm.dest( toRank, X0Z1) );

      /*
        MPI_Isend(destAddr, xferFields*dy, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]);
      */
      ++emsg;
    }

    if (rowNotMin && colNotMin && planeNotMin) {
      // corner at domain logical coord (0, 0, 0)
      int toRank = myRank - domain.tp()*domain.tp() - domain.tp() - 1;
      Real_t *comBuf = &comm.commDataSend()[pmsg * maxPlaneComm +
                                            emsg * maxEdgeComm +
                                            cmsg * CACHE_COHERENCE_PAD_REAL];
      for (Index_t fi=0; fi<xferFields; ++fi) {
        comBuf[fi] = (domain.*fieldData[fi])(0);
      }

      comm.sendRequest[pmsg + emsg + cmsg] =
        dash::copy_async(
          comBuf, comBuf+xferFields,
          comm.dest( toRank, X1Y1Z1) );
      /*
        MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]);
      */
      ++cmsg;
    }
    if (rowNotMin && colNotMin && planeNotMax && doSend) {
      // corner at domain logical coord (0, 0, 1)
      int toRank = myRank + domain.tp()*domain.tp() - domain.tp() - 1;
      Real_t *comBuf = &comm.commDataSend()[pmsg * maxPlaneComm +
                                            emsg * maxEdgeComm +
                                            cmsg * CACHE_COHERENCE_PAD_REAL];
      Index_t idx = dx*dy*(dz - 1);
      for (Index_t fi=0; fi<xferFields; ++fi) {
        comBuf[fi] = (domain.*fieldData[fi])(idx);
      }

      comm.sendRequest[pmsg + emsg + cmsg] =
        dash::copy_async(
          comBuf, comBuf+xferFields,
          comm.dest( toRank, X1Y1Z0 ) );
      /*
        MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]);
      */
      ++cmsg;
    }
    if (rowNotMin && colNotMax && planeNotMin) {
      // corner at domain logical coord (1, 0, 0)
      int toRank = myRank - domain.tp()*domain.tp() - domain.tp() + 1;
      Real_t *comBuf = &comm.commDataSend()[pmsg * maxPlaneComm +
                                            emsg * maxEdgeComm +
                                            cmsg * CACHE_COHERENCE_PAD_REAL];
      Index_t idx = dx - 1;
      for (Index_t fi=0; fi<xferFields; ++fi) {
        comBuf[fi] = (domain.*fieldData[fi])(idx);
      }

      comm.sendRequest[pmsg + emsg + cmsg] =
        dash::copy_async(
          comBuf, comBuf+xferFields,
          comm.dest( toRank, X0Y1Z1) );
      /*
        MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]);
      */
      ++cmsg;
    }
    if (rowNotMin && colNotMax && planeNotMax && doSend) {
      // corner at domain logical coord (1, 0, 1)
      int toRank = myRank + domain.tp()*domain.tp() - domain.tp() + 1;
      Real_t *comBuf = &comm.commDataSend()[pmsg * maxPlaneComm +
                                            emsg * maxEdgeComm +
                                            cmsg * CACHE_COHERENCE_PAD_REAL];
      Index_t idx = dx*dy*(dz - 1) + (dx - 1);
      for (Index_t fi=0; fi<xferFields; ++fi) {
        comBuf[fi] = (domain.*fieldData[fi])(idx);
      }

      comm.sendRequest[pmsg + emsg + cmsg] =
        dash::copy_async(
          comBuf, comBuf+xferFields,
          comm.dest( toRank, X0Y1Z0 ) );
      /*
        MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]);
      */
      ++cmsg;
    }
    if (rowNotMax && colNotMin && planeNotMin) {
      // corner at domain logical coord (0, 1, 0)
      int toRank = myRank - domain.tp()*domain.tp() + domain.tp() - 1;
      Real_t *comBuf = &comm.commDataSend()[pmsg * maxPlaneComm +
                                            emsg * maxEdgeComm +
                                            cmsg * CACHE_COHERENCE_PAD_REAL];
      Index_t idx = dx*(dy - 1);
      for (Index_t fi=0; fi<xferFields; ++fi) {
        comBuf[fi] = (domain.*fieldData[fi])(idx);
      }

      comm.sendRequest[pmsg + emsg + cmsg] =
        dash::copy_async(
          comBuf, comBuf+xferFields,
          comm.dest( toRank, X1Y0Z1) );
      /*
        MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]);
      */
      ++cmsg;
    }
    if (rowNotMax && colNotMin && planeNotMax && doSend) {
      // corner at domain logical coord (0, 1, 1)
      int toRank = myRank + domain.tp()*domain.tp() + domain.tp() - 1;
      Real_t *comBuf = &comm.commDataSend()[pmsg * maxPlaneComm +
                                            emsg * maxEdgeComm +
                                            cmsg * CACHE_COHERENCE_PAD_REAL];
      Index_t idx = dx*dy*(dz - 1) + dx*(dy - 1);
      for (Index_t fi=0; fi<xferFields; ++fi) {
        comBuf[fi] = (domain.*fieldData[fi])(idx);
      }

      comm.sendRequest[pmsg + emsg + cmsg] =
        dash::copy_async(
          comBuf, comBuf+xferFields,
          comm.dest( toRank, X1Y0Z0 ) );
      /*
        MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]);
      */
      ++cmsg;
    }
    if (rowNotMax && colNotMax && planeNotMin) {
      // corner at domain logical coord (1, 1, 0)
      int toRank = myRank - domain.tp()*domain.tp() + domain.tp() + 1;
      Real_t *comBuf = &comm.commDataSend()[pmsg * maxPlaneComm +
                                            emsg * maxEdgeComm +
                                            cmsg * CACHE_COHERENCE_PAD_REAL];
      Index_t idx = dx*dy - 1;
      for (Index_t fi=0; fi<xferFields; ++fi) {
        comBuf[fi] = (domain.*fieldData[fi])(idx);
      }

      comm.sendRequest[pmsg + emsg + cmsg] =
        dash::copy_async(
          comBuf, comBuf+xferFields,
          comm.dest( toRank, X0Y0Z1 ) );
      /*
        MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
        MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]);
      */
      ++cmsg;
    }
    if (rowNotMax && colNotMax && planeNotMax && doSend) {
      // corner at domain logical coord (1, 1, 1)
      int toRank = myRank + domain.tp()*domain.tp() + domain.tp() + 1;
      Real_t *comBuf = &comm.commDataSend()[pmsg * maxPlaneComm +
                                            emsg * maxEdgeComm +
                                            cmsg * CACHE_COHERENCE_PAD_REAL];
      Index_t idx = dx*dy*dz - 1;
      for (Index_t fi=0; fi<xferFields; ++fi) {
        comBuf[fi] = (domain.*fieldData[fi])(idx);
      }

      comm.sendRequest[pmsg + emsg + cmsg] =
        dash::copy_async(
          comBuf, comBuf+xferFields,
          comm.dest( toRank, X0Y0Z0 ) );
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

  int myRank;

  Index_t dx = domain.sizeX() + 1;
  Index_t dy = domain.sizeY() + 1;
  Index_t dz = domain.sizeZ() + 1;

  Real_t *srcAddr;
  Index_t rowNotMin, rowNotMax, colNotMin, colNotMax, planeNotMin, planeNotMax;

  // assume communication to 6 neighbors by default
  rowNotMin = rowNotMax = colNotMin = colNotMax = planeNotMin = planeNotMax = true;

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

    if (planeNotMin) {
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
    }
    if (planeNotMax) {
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
    }
  }

  if (rowNotMin | rowNotMax) {
    /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
    Index_t opCount = dx * dz;

    if (rowNotMin) {
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
    }
    if (rowNotMax) {
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
    }
  }
  if (colNotMin | colNotMax) {
    /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
    Index_t opCount = dy * dz;

    if (colNotMin) {
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
    }
    if (colNotMax) {
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
    }
  }

  if (rowNotMin & colNotMin) {
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
  }

  if (rowNotMin & planeNotMin) {
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
  }

  if (colNotMin & planeNotMin) {
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
  }

  if (rowNotMax & colNotMax) {
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
  }

  if (rowNotMax & planeNotMax) {
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
  }

  if (colNotMax & planeNotMax) {
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
  }

  if (rowNotMax & colNotMin) {
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
  }

  if (rowNotMin & planeNotMax) {
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
  }

  if (colNotMin & planeNotMax) {
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
  }

  if (rowNotMin & colNotMax) {
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
  }

  if (rowNotMax & planeNotMin) {
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
  }

  if (colNotMax & planeNotMin) {
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
  }

  if (rowNotMin & colNotMin & planeNotMin) {
    /* corner at domain logical coord (0, 0, 0) */
    Real_t *comBuf = &comm.commDataRecv()[comm.offset(X0Y0Z0)];
    DBGSYNC(xferFields, 1, X0Y0Z0);
    //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      (domain.*fieldData[fi])(0) += comBuf[fi];
    }
  }
  if (rowNotMin & colNotMin & planeNotMax) {
    /* corner at domain logical coord (0, 0, 1) */
    Real_t *comBuf = &comm.commDataRecv()[comm.offset(X0Y0Z1)];
    DBGSYNC(xferFields, 1, X0Y0Z1);
    Index_t idx = dx*dy*(dz - 1);
    //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      (domain.*fieldData[fi])(idx) += comBuf[fi];
    }
  }
  if (rowNotMin & colNotMax & planeNotMin) {
    /* corner at domain logical coord (1, 0, 0) */
    Real_t *comBuf = &comm.commDataRecv()[comm.offset(X1Y0Z0)];
    DBGSYNC(xferFields, 1, X1Y0Z0);
    Index_t idx = dx - 1;
    //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      (domain.*fieldData[fi])(idx) += comBuf[fi];
    }
  }
  if (rowNotMin & colNotMax & planeNotMax) {
    /* corner at domain logical coord (1, 0, 1) */
    Real_t *comBuf = &comm.commDataRecv()[comm.offset(X1Y0Z1)];
    DBGSYNC(xferFields, 1, X1Y0Z1);
    Index_t idx = dx*dy*(dz - 1) + (dx - 1);
    //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      (domain.*fieldData[fi])(idx) += comBuf[fi];
    }
  }
  if (rowNotMax & colNotMin & planeNotMin) {
    /* corner at domain logical coord (0, 1, 0) */
    Real_t *comBuf = &comm.commDataRecv()[comm.offset(X0Y1Z0)];
    DBGSYNC(xferFields, 1, X0Y1Z0);
    Index_t idx = dx*(dy - 1);
    //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      (domain.*fieldData[fi])(idx) += comBuf[fi];
    }
  }
  if (rowNotMax & colNotMin & planeNotMax) {
    /* corner at domain logical coord (0, 1, 1) */
    Real_t *comBuf = &comm.commDataRecv()[comm.offset(X0Y1Z1)];
    DBGSYNC(xferFields, 1, X0Y1Z1);
    Index_t idx = dx*dy*(dz - 1) + dx*(dy - 1);
    //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      (domain.*fieldData[fi])(idx) += comBuf[fi];
    }
  }
  if (rowNotMax & colNotMax & planeNotMin) {
    /* corner at domain logical coord (1, 1, 0) */
    Real_t *comBuf = &comm.commDataRecv()[comm.offset(X1Y1Z0)];
    DBGSYNC(xferFields, 1, X1Y1Z0);
    Index_t idx = dx*dy - 1;
    //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      (domain.*fieldData[fi])(idx) += comBuf[fi];
    }
  }
  if (rowNotMax & colNotMax & planeNotMax) {
    /* corner at domain logical coord (1, 1, 1) */
    Real_t *comBuf = &comm.commDataRecv()[comm.offset(X1Y1Z1)];
    DBGSYNC(xferFields, 1, X1Y1Z1);
    Index_t idx = dx*dy*dz - 1;
    //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      (domain.*fieldData[fi])(idx) += comBuf[fi];
    }
  }
}

void DASHCommSyncPosVel(Domain& domain, DASHComm& comm)
{
  if (domain.numRanks() == 1)
    return;

  int myRank;
  bool doRecv = false;
  Index_t xferFields = 6; /* x, y, z, xd, yd, zd */
  Domain_member fieldData[6];
  Index_t maxPlaneComm = xferFields * domain.maxPlaneSize();
  Index_t maxEdgeComm  = xferFields * domain.maxEdgeSize();
  Index_t dx = domain.sizeX() + 1;
  Index_t dy = domain.sizeY() + 1;
  Index_t dz = domain.sizeZ() + 1;
  //MPI_Status status;
  Real_t *srcAddr;
  bool rowNotMin, rowNotMax, colNotMin, colNotMax, planeNotMin, planeNotMax;
  // assume communication to 6 neighbors by default
  rowNotMin = rowNotMax = colNotMin = colNotMax = planeNotMin = planeNotMax = true;

  /* assume communication to 6 neighbors by default */
  if( domain.rowLoc()   == 0 )               { rowNotMin   = false; }
  if( domain.rowLoc()   == (domain.tp()-1) ) { rowNotMax   = false; }
  if( domain.colLoc()   == 0 )               { colNotMin   = false; }
  if( domain.colLoc()   == (domain.tp()-1) ) { colNotMax   = false; }
  if( domain.planeLoc() == 0 )               { planeNotMin = false; }
  if( domain.planeLoc() == (domain.tp()-1) ) { planeNotMax = false; }

  fieldData[0] = &Domain::x;
  fieldData[1] = &Domain::y;
  fieldData[2] = &Domain::z;
  fieldData[3] = &Domain::xd;
  fieldData[4] = &Domain::yd;
  fieldData[5] = &Domain::zd;

  myRank = dash::myid();
  if (planeNotMin | planeNotMax) {
    /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
    Index_t opCount = dx * dy;

    if (planeNotMin && doRecv) {
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
    }
    if (planeNotMax) {
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
    }
  }
  if (rowNotMin | rowNotMax) {
    /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
    Index_t opCount = dx * dz;

    if (rowNotMin && doRecv) {
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
    }
    if (rowNotMax) {
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
    }
  }

  if (colNotMin | colNotMax) {
    /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
    Index_t opCount = dy * dz;

    if (colNotMin && doRecv) {
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
    }
    if (colNotMax) {
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
    }
  }
  if (rowNotMin && colNotMin && doRecv) {
    srcAddr = &comm.commDataRecv()[comm.offset(X0Y0)];
    //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      Domain_member dest = fieldData[fi];
      for (Index_t i=0; i<dz; ++i) {
        (domain.*dest)(i*dx*dy) = srcAddr[i];
      }
      srcAddr += dz;
    }
  }

  if (rowNotMin && planeNotMin && doRecv) {
    srcAddr = &comm.commDataRecv()[comm.offset(X0Y0)];
    //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      Domain_member dest = fieldData[fi];
      for (Index_t i=0; i<dx; ++i) {
        (domain.*dest)(i) = srcAddr[i];
      }
      srcAddr += dx;
    }
  }

  if (colNotMin && planeNotMin && doRecv) {
    srcAddr = &comm.commDataRecv()[comm.offset(X0Z0)];
    //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      Domain_member dest = fieldData[fi];
      for (Index_t i=0; i<dy; ++i) {
        (domain.*dest)(i*dx) = srcAddr[i];
      }
      srcAddr += dy;
    }
  }

  if (rowNotMax && colNotMax) {
    srcAddr = &comm.commDataRecv()[comm.offset(X1Y1)];
    //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      Domain_member dest = fieldData[fi];
      for (Index_t i=0; i<dz; ++i) {
        (domain.*dest)(dx*dy - 1 + i*dx*dy) = srcAddr[i];
      }
      srcAddr += dz;
    }
  }

  if (rowNotMax && planeNotMax) {
    srcAddr = &comm.commDataRecv()[comm.offset(Y1Z1)];
    //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      Domain_member dest = fieldData[fi];
      for (Index_t i=0; i<dx; ++i) {
        (domain.*dest)(dx*(dy-1) + dx*dy*(dz-1) + i) = srcAddr[i];
      }
      srcAddr += dx;
    }
  }

  if (colNotMax && planeNotMax) {
    srcAddr = &comm.commDataRecv()[comm.offset(X1Z1)];
    //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      Domain_member dest = fieldData[fi];
      for (Index_t i=0; i<dy; ++i) {
        (domain.*dest)(dx*dy*(dz-1) + dx - 1 + i*dx) = srcAddr[i];
      }
      srcAddr += dy;
    }
  }

  if (rowNotMax && colNotMin) {
    srcAddr = &comm.commDataRecv()[comm.offset(X0Y1)];
    //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      Domain_member dest = fieldData[fi];
      for (Index_t i=0; i<dz; ++i) {
        (domain.*dest)(dx*(dy-1) + i*dx*dy) = srcAddr[i];
      }
      srcAddr += dz;
    }
  }

  if (rowNotMin && planeNotMax) {
    srcAddr = &comm.commDataRecv()[comm.offset(Y0Z1)];
    //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      Domain_member dest = fieldData[fi];
      for (Index_t i=0; i<dx; ++i) {
        (domain.*dest)(dx*dy*(dz-1) + i) = srcAddr[i];
      }
      srcAddr += dx;
    }
  }

  if (colNotMin && planeNotMax) {
    srcAddr = &comm.commDataRecv()[comm.offset(X0Z1)];
    //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      Domain_member dest = fieldData[fi];
      for (Index_t i=0; i<dy; ++i) {
        (domain.*dest)(dx*dy*(dz-1) + i*dx) = srcAddr[i];
      }
      srcAddr += dy;
    }
  }

  if (rowNotMin && colNotMax && doRecv) {
    srcAddr = &comm.commDataRecv()[comm.offset(X1Y0)];
    //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      Domain_member dest = fieldData[fi];
      for (Index_t i=0; i<dz; ++i) {
        (domain.*dest)(dx - 1 + i*dx*dy) = srcAddr[i];
      }
      srcAddr += dz;
    }
  }

  if (rowNotMax && planeNotMin && doRecv) {
    srcAddr = &comm.commDataRecv()[comm.offset(Y1Z0)];
    //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      Domain_member dest = fieldData[fi];
      for (Index_t i=0; i<dx; ++i) {
        (domain.*dest)(dx*(dy - 1) + i) = srcAddr[i];
      }
      srcAddr += dx;
    }
  }

  if (colNotMax && planeNotMin && doRecv) {
    srcAddr = &comm.commDataRecv()[comm.offset(X1Z0)];
    //MPI_Wait(&comm.recvRequest[pmsg+emsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      Domain_member dest = fieldData[fi];
      for (Index_t i=0; i<dy; ++i) {
        (domain.*dest)(dx - 1 + i*dx) = srcAddr[i];
      }
      srcAddr += dy;
    }
  }

  if (rowNotMin && colNotMin && planeNotMin && doRecv) {
    /* corner at domain logical coord (0, 0, 0) */
    Real_t *comBuf = &comm.commDataRecv()[comm.offset(X0Y0Z0)];
    //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      (domain.*fieldData[fi])(0) = comBuf[fi];
    }
  }
  if (rowNotMin && colNotMin && planeNotMax) {
    /* corner at domain logical coord (0, 0, 1) */
    Real_t *comBuf = &comm.commDataRecv()[comm.offset(X0Y0Z1)];
    Index_t idx = dx*dy*(dz - 1);
    //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      (domain.*fieldData[fi])(idx) = comBuf[fi];
    }
  }
  if (rowNotMin && colNotMax && planeNotMin && doRecv) {
    /* corner at domain logical coord (1, 0, 0) */
    Real_t *comBuf = &comm.commDataRecv()[comm.offset(X1Y0Z0)];
    Index_t idx = dx - 1;
    //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      (domain.*fieldData[fi])(idx) = comBuf[fi];
    }
  }
  if (rowNotMin && colNotMax && planeNotMax) {
    /* corner at domain logical coord (1, 0, 1) */
    Real_t *comBuf = &comm.commDataRecv()[comm.offset(X1Y0Z1)];
    Index_t idx = dx*dy*(dz - 1) + (dx - 1);
    //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      (domain.*fieldData[fi])(idx) = comBuf[fi];
    }
  }
  if (rowNotMax && colNotMin && planeNotMin && doRecv) {
    /* corner at domain logical coord (0, 1, 0) */
    Real_t *comBuf = &comm.commDataRecv()[comm.offset(X0Y1Z0)];
    Index_t idx = dx*(dy - 1);
    //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      (domain.*fieldData[fi])(idx) = comBuf[fi];
    }
  }
  if (rowNotMax && colNotMin && planeNotMax) {
    /* corner at domain logical coord (0, 1, 1) */
    Real_t *comBuf = &comm.commDataRecv()[comm.offset(X0Y1Z1)];
    Index_t idx = dx*dy*(dz - 1) + dx*(dy - 1);
    //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      (domain.*fieldData[fi])(idx) = comBuf[fi];
    }
  }
  if (rowNotMax && colNotMax && planeNotMin && doRecv) {
    /* corner at domain logical coord (1, 1, 0) */
    Real_t *comBuf = &comm.commDataRecv()[comm.offset(X1Y1Z0)];
    Index_t idx = dx*dy - 1;
    //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      (domain.*fieldData[fi])(idx) = comBuf[fi];
    }
  }
  if (rowNotMax && colNotMax && planeNotMax) {
    /* corner at domain logical coord (1, 1, 1) */
    Real_t *comBuf = &comm.commDataRecv()[comm.offset(X1Y1Z1)];
    Index_t idx = dx*dy*dz - 1;
    //MPI_Wait(&comm.recvRequest[pmsg+emsg+cmsg], &status);
    for (Index_t fi=0; fi<xferFields; ++fi) {
      (domain.*fieldData[fi])(idx) = comBuf[fi];
    }
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
  Real_t *srcAddr;
  bool rowNotMin, rowNotMax, colNotMin, colNotMax, planeNotMin, planeNotMax;
  // assume communication to 6 neighbors by default
  rowNotMin = rowNotMax = colNotMin = colNotMax = planeNotMin = planeNotMax = true;

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
      /* contiguous memory */
      srcAddr = &comm.commDataRecv()[comm.offset(Z0)];
      //MPI_Wait(&comm.recvRequest[pmsg], &status);
      for (Index_t fi=0; fi<xferFields; ++fi) {
        Domain_member dest = fieldData[fi];
        for (Index_t i=0; i<opCount; ++i) {
          (domain.*dest)(fieldOffset[fi] + i) = srcAddr[i];
        }
        srcAddr += opCount;
        fieldOffset[fi] += opCount;
      }
    }
    if (planeNotMax) {
      /* contiguous memory */
      srcAddr = &comm.commDataRecv()[comm.offset(Z1)];
      //MPI_Wait(&comm.recvRequest[pmsg], &status);
      for (Index_t fi=0; fi<xferFields; ++fi) {
        Domain_member dest = fieldData[fi];
        for (Index_t i=0; i<opCount; ++i) {
          (domain.*dest)(fieldOffset[fi] + i) = srcAddr[i];
        }
        srcAddr += opCount;
        fieldOffset[fi] += opCount;
      }
    }
  }

  if (rowNotMin | rowNotMax) {
    /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
    Index_t opCount = dx * dz;

    if (rowNotMin) {
      /* contiguous memory */
      srcAddr = &comm.commDataRecv()[comm.offset(Y0)];
      //MPI_Wait(&comm.recvRequest[pmsg], &status);
      for (Index_t fi=0; fi<xferFields; ++fi) {
        Domain_member dest = fieldData[fi];
        for (Index_t i=0; i<opCount; ++i) {
          (domain.*dest)(fieldOffset[fi] + i) = srcAddr[i];
        }
        srcAddr += opCount;
        fieldOffset[fi] += opCount;
      }
    }
    if (rowNotMax) {
      /* contiguous memory */
      srcAddr = &comm.commDataRecv()[comm.offset(Y1)];
      //MPI_Wait(&comm.recvRequest[pmsg], &status);
      for (Index_t fi=0; fi<xferFields; ++fi) {
        Domain_member dest = fieldData[fi];
        for (Index_t i=0; i<opCount; ++i) {
          (domain.*dest)(fieldOffset[fi] + i) = srcAddr[i];
        }
        srcAddr += opCount;
        fieldOffset[fi] += opCount;
      }
    }
  }
  if (colNotMin | colNotMax) {
    /* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
    Index_t opCount = dy * dz;

    if (colNotMin) {
      /* contiguous memory */
      srcAddr = &comm.commDataRecv()[comm.offset(X0)];
      //MPI_Wait(&comm.recvRequest[pmsg], &status);
      for (Index_t fi=0; fi<xferFields; ++fi) {
        Domain_member dest = fieldData[fi];
        for (Index_t i=0; i<opCount; ++i) {
          (domain.*dest)(fieldOffset[fi] + i) = srcAddr[i];
        }
        srcAddr += opCount;
        fieldOffset[fi] += opCount;
      }
    }
    if (colNotMax) {
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
    }
  }
}
