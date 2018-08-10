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
 * CODE.H: define various global things for CODE.C.
 */

#ifndef _CODE_H_
#define _CODE_H_

#include <dash/Array.h>
#include <dash/Mutex.h>
#include <dash/Shared.h>
#include <type_traits>
#include "defs.h"
#include "shared_array.h"
#include "stdinc.h"
#include <cassert>

#define PAD_SIZE (PAGE_SIZE / (sizeof(long)))
/* Workaround since dash::Shared does not support arrays as value types */
typedef struct {
  real x0y0, x0y1, x0y2;
  real x1y0, x1y1, x1y2;
  real x2y0, x2y1, x2y2;
} sh_mat;

typedef struct {
  real x, y, z;
} sh_vec;

/* Defined by the input file */
// string headline;   /* message describing calculation */
// string infile;     /* file name for snapshot input */
// string outfile;    /* file name for snapshot output */
extern size_t nbody;   /* number of bodies in system */
extern real   dtime;   /* timestep for leapfrog integrator */
extern real   dtout;   /* time between data outputs */
extern real   tstop;   /* time to stop calculation */
extern real   fcells;  /* ratio of cells/leaves allocated */
extern real   fleaves; /* ratio of leaves/bodies allocated */
extern real   tol;     /* accuracy parameter: 0.0 => exact */
extern real   tolsq;   /* square of previous */
extern real   eps;     /* potential softening parameter */
extern real   epssq;   /* square of previous */
extern real   dthf;    /* half time step */
// dash::Shared<long> NPROC;    /* Number of Processors */

extern long   maxcell;   /* max number of cells allocated */
extern long   maxleaf;   /* max number of leaves allocated */
extern size_t maxmybody; /* max no. of bodies allocated per processor */
extern size_t maxmycell; /* max num. of cells to be allocated */
extern size_t maxmyleaf; /* max num. of leaves to be allocated */

extern dash::Array<body> bodytab; /* array size is exactly nbody bodies */
extern dash::Array<cell> celltab;
extern dash::Array<leaf> leaftab;

extern std::vector<dash::Mutex> CellLock;

// struct GlobalMemory  {  /* all this info is for the whole system */

/* total number of body/cell interactions  */
// extern dash::Shared<long> n2bcalc;
/* total number of body/body interactions  */
// extern dash::Shared<long> nbccalc;
/* number of self interactions */
// extern dash::Shared<long> selfint;

extern dash::Shared<real> mtot; /* total mass of N-body system             */
// dash::Shared<real> etot[3];      /* binding, kinetic, potential energy */
// extern dash::Shared<sh_mat> keten; /* kinetic energy tensor */
// extern dash::Shared<sh_mat> peten; /* potential energy tensor */
// dash::Shared<vector> cmphase[2]; /* center of mass coordinates and velocity
// */
// extern dash::Shared<sh_vec>  amvec;  /* angular momentum vector */
extern dash::Shared<cellptr> G_root; /* root of the whole tree */
/* lower-left corner of coordinate box */
extern vector g_rmin;
/* temporary lower-left corner of the box  */
extern vector g_min;
/* temporary upper right corner of the box */
extern vector g_max;
/* side-length of integer coordinate box   */
extern real g_rsize;

extern struct local_memory Local;

// pthread_barrier_t (Barrier);

/* barrier at the beginning of stepsystem  */
#if 0
  pthread_mutex_t (CountLock); /* Lock on the shared variables            */
  pthread_mutex_t (NcellLock); /* Lock on the counter of array of cells for loadtree */
  pthread_mutex_t (NleafLock);/* Lock on the counter of array of leaves for loadtree */
  pthread_mutex_t (io_lock);
#endif
/*
dash::Shared<unsigned long> createstart,createend,computestart,computeend;
dash::Shared<unsigned long> trackstart, trackend, tracktime;
dash::Shared<unsigned long> partitionstart, partitionend, partitiontime;
dash::Shared<unsigned long> treebuildstart, treebuildend, treebuildtime;
dash::Shared<unsigned long> forcecalcstart, forcecalcend, forcecalctime;
*/
// long current_id;
// volatile long k; /*for memory allocation in code.C */
//};
// global struct GlobalMemory *Global;

/* This structure is needed because under the sproc model there is no
 * per processor private address space.
 */
struct local_memory {
  /* Use padding so that each processor's variables are on their own page */
  // long pad_begin[PAD_SIZE];

  real tnow;  /* current value of simulation time */
  real tout;  /* time next output is due */
  long nstep; /* number of integration steps so far */

  long workMin, workMax; /* interval of cost to be treated by a proc */

  vector min = {0}, max = {0}; /* min and max of coordinates for each Proc. */

  /* num. of cells used for this proc in ctab */
  size_t mynumcell;
  /* num. of leaves used for this proc in ctab */
  size_t mynumleaf;
  /* num bodies allocated to the processor */
  size_t mynbody;
  /* array of bodies allocated / processor */
  std::vector<bodyptr> mybodytab;
  /* num cells allocated to the processor */
  long myncell;
  /* array of cellptrs allocated to the processor */
  std::vector<cellptr> mycelltab;
  /* number of leaves allocated to the processor */
  long mynleaf;
  /* array of leafptrs allocated to the processor */
  std::vector<leafptr> myleaftab;

  // decltype(celltab)::local_type ctab;  /* array of cells used for the tree.
  // */
  // decltype(leaftab)::local_type ltab;  /* array of cells used for the tree.
  // */

  long    myn2bcalc;  /* body-body force calculations for each processor */
  long    mynbccalc;  /* body-cell force calculations for each processor */
  long    myselfint;  /* count self-interactions for each processor */
  long    myn2bterm;  /* count body-body terms for a body */
  long    mynbcterm;  /* count body-cell terms for a body */
  bool    skipself;   /* true if self-interaction skipped OK */
  bodyptr pskip;      /* body to skip in force evaluation */
  vector  pos0 = {0}; /* point at which to evaluate field */
  real    phi0 = {0}; /* computed potential at pos0 */
  vector  acc0 = {0}; /* computed acceleration at pos0 */
  vector  dr   = {0}; /* data to be shared */
  real    drsq;       /* between gravsub and subdivp */
  nodeptr pmem;       /* remember particle data */

  nodeptr Current_Root;  // Cell pointer
  long    Root_Coords[NDIM];

  real    mymtot;        /* total mass of N-body system */
  real    myetot[3];     /* binding, kinetic, potential energy */
  matrix  myketen;       /* kinetic energy tensor */
  matrix  mypeten;       /* potential energy tensor */
  vector  mycmphase[2];  /* center of mass coordinates */
  vector  myamvec = {0}; /* angular momentum vector */
  uint8_t pad[PAGE_SIZE];
};

void SlaveStart(void);
void stepsystem(long ProcessId);
void ComputeForces(long ProcessId);
void Help(void);
void ANLinit(void);
void init_root(void);
void tab_init(void);
void startrun(void);
void testdata(void);
void pickshell(real* vec, real rad);
void find_my_initial_bodies(
    dash::Array<body>& btab, long nbody, long ProcessId);
void find_my_bodies(
    nodeptr mycell, long work, long direction, long ProcessId);
void Housekeep(long ProcessId);
void setbound(void);
long Log_base_2(long number);
long intpow(long i, long j);

#ifndef NDEBUG
#define ASSERT(expr) (assert((expr)))
#else
#define ASSERT(expr) (void)(expr)
#endif

#endif
