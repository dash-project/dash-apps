/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/

/*
 * GRAV.C:
 */

#include <pthread.h>

#include <sys/time.h>

#include <unistd.h>

#include <stdlib.h>

#include <malloc.h>

#include "code.h"
#include "defs.h"
#include "grav.h"

#include "stdinc.h"

/*
 * HACKGRAV: evaluate grav field at a given particle.
 */

void hackgrav(bodyptr p, long ProcessId)
{
  // Skip the body itself as a body cannot apply some forces to itself
  Local.pskip = p;
  body p_val  = *p;
  SETV(Local.pos0, p_val.pos);
  Local.phi0 = 0.0;
  CLRV(Local.acc0);
  Local.myn2bterm = 0;
  Local.mynbcterm = 0;
  Local.skipself  = FALSE;
  hackwalk(ProcessId);
  p_val.phi = Local.phi0;
  SETV(p_val.acc, Local.acc0);
#ifdef QUADPOLE
  p_val.cost = Local.myn2bterm + NDIM * Local.mynbcterm;
#else
  p_val.cost = Local.myn2bterm + Local.mynbcterm;
#endif
  *p = p_val;
}

/*
 * GRAVSUB: compute a single body-body or body-cell longeraction.
 */

void gravsub(nodeptr p, long ProcessId)
{
  real drabs, phii, mor3;
  vector ai;

  auto const p_val = static_cast<node>(*p);

  if (p != Local.pmem) {
    SUBV(Local.dr, p_val.pos, Local.pos0);
    DOTVP(Local.drsq, Local.dr, Local.dr);
  }

  Local.drsq += epssq.get();
  drabs = sqrt((double)Local.drsq);
  phii  = p_val.mass / drabs;
  Local.phi0 -= phii;
  mor3 = phii / Local.drsq;
  MULVS(ai, Local.dr, mor3);
  ADDV(Local.acc0, Local.acc0, ai);
  if (p_val.type != BODY) { /* a body-cell/leaf interaction? */
    Local.mynbcterm++;
#ifdef QUADPOLE
    dr5inv = 1.0 / (Local.drsq * Local.drsq * drabs);
    MULMV(quaddr, Quad(p), Local.dr);
    DOTVP(drquaddr, Local.dr, quaddr);
    phiquad = -0.5 * dr5inv * drquaddr;
    Local.phi0 += phiquad;
    phiquad = 5.0 * phiquad / Local.drsq;
    MULVS(ai, Local.dr, phiquad);
    SUBV(Local.acc0, Local.acc0, ai);
    MULVS(quaddr, quaddr, dr5inv);
    SUBV(Local.acc0, Local.acc0, quaddr);
#endif
  }
  else { /* a body-body interaction  */
    Local.myn2bterm++;
  }
}

/*
 * HACKWALK: walk the tree opening cells too close to a given point.
 */

void hackwalk(long ProcessId)
{
  walksub(G_root.get().get(),
          rsize.get().get() * rsize.get().get(), ProcessId);
}

/*
 * WALKSUB: recursive routine to do hackwalk operation.
 */

void walksub(nodeptr n, real dsq, long ProcessId)
{
  leafptr l;
  bodyptr p;
  long i;

  if (subdivp(n, dsq, ProcessId)) {
    if (Type(*n) == CELL) {
      auto const n_val = static_cast<cell>(*NODE_AS_CELL(n));

      for (auto nn = n_val.subp; nn < n_val.subp + NSUB; nn++) {
        if (*nn != nodeptr(nullptr)) {
          walksub(*nn, dsq / 4.0, ProcessId);
        }
      }
    }
    else {
      l = NODE_AS_LEAF(n);
      auto const l_val = static_cast<leaf>(*l);
      for (i = 0; i < l_val.num_bodies; i++) {
        p = l_val.bodyp[i];
        if (p != Local.pskip) {
          gravsub(p, ProcessId);
        }
        else {
          Local.skipself = TRUE;
        }
      }
    }
  }
  else {
    gravsub(n, ProcessId);
  }
}

/*
 * SUBDIVP: decide if a node should be opened.
 * Side effects: sets  pmem,dr, and drsq.
 */

bool subdivp(nodeptr const p, real dsq, long ProcessId)
{
  ASSERT(p);
  auto const p_val = static_cast<node>(*p);
  SUBV(Local.dr, p_val.pos, Local.pos0);
  DOTVP(Local.drsq, Local.dr, Local.dr);
  Local.pmem = p;
  return (tolsq.get() * Local.drsq < dsq);
}
