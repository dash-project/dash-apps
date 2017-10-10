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

#include "stdinc.h"
#include "defs.h"
#include <dash/Shared.h>
#include <dash/Array.h>
#include <dash/Mutex.h>
#include <type_traits>

#define PAD_SIZE (PAGE_SIZE / (sizeof(long)))

/* Defined by the input file */
//string headline;   /* message describing calculation */
 //string infile;     /* file name for snapshot input */
 //string outfile;    /* file name for snapshot output */
extern long nbody;    /* number of bodies in system */
extern dash::Shared<real> dtime;    /* timestep for leapfrog integrator */
extern dash::Shared<real> dtout;    /* time between data outputs */
extern dash::Shared<real> tstop;    /* time to stop calculation */
extern dash::Shared<real> fcells;     /* ratio of cells/leaves allocated */
extern dash::Shared<real> fleaves;    /* ratio of leaves/bodies allocated */
extern dash::Shared<real> tol;    /* accuracy parameter: 0.0 => exact */
extern dash::Shared<real> tolsq;    /* square of previous */
extern dash::Shared<real> eps;    /* potential softening parameter */
extern dash::Shared<real> epssq;    /* square of previous */
extern dash::Shared<real> dthf;     /* half time step */
// dash::Shared<long> NPROC;    /* Number of Processors */

extern dash::Shared<long> maxcell; /* max number of cells allocated */
extern dash::Shared<long> maxleaf; /* max number of leaves allocated */
extern long maxmybody;             /* max no. of bodies allocated per processor */
extern long maxmycell;             /* max num. of cells to be allocated */
extern long maxmyleaf;             /* max num. of leaves to be allocated */

extern dash::Array<body> bodytab; /* array size is exactly nbody bodies */
extern dash::Array<cell> celltab;
extern dash::Array<leaf> leaftab;

extern std::vector<dash::Mutex> CellLock;
extern dash::Mutex CountLock;

typedef struct {
  real x0y0, x0y1, x0y2;
  real x1y0, x1y1, x1y2;
  real x2y0, x2y1, x2y2;
} sh_mat;

typedef struct {
  real x, y, z;
} sh_vec;


//struct GlobalMemory  {  /* all this info is for the whole system */
  extern dash::Shared<long> n2bcalc; /* total number of body/cell interactions  */
  extern dash::Shared<long> nbccalc; /* total number of body/body interactions  */
  extern dash::Shared<long> selfint; /* number of self interactions             */
  extern dash::Shared<real> mtot;    /* total mass of N-body system             */
  // dash::Shared<real> etot[3];      /* binding, kinetic, potential energy */
  extern dash::Shared<sh_mat> keten; /* kinetic energy tensor                   */
  extern dash::Shared<sh_mat> peten; /* potential energy tensor                 */
  // dash::Shared<vector> cmphase[2]; /* center of mass coordinates and velocity
  // */
  extern dash::Shared<sh_vec> amvec;   /* angular momentum vector                 */
  extern dash::Shared<cellptr> G_root; /* root of the whole tree                  */
  extern dash::Shared<sh_vec> rmin;    /* lower-left corner of coordinate box     */
  extern dash::Shared<sh_vec> min;     /* temporary lower-left corner of the box  */
  extern dash::Shared<sh_vec> max;     /* temporary upper right corner of the box */
  extern dash::Shared<real> rsize;     /* side-length of integer coordinate box   */

//pthread_barrier_t (Barrier);

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
  //long current_id;
  //volatile long k; /*for memory allocation in code.C */
//};
//global struct GlobalMemory *Global;

/* This structure is needed because under the sproc model there is no
 * per processor private address space.
 */
struct local_memory {
   /* Use padding so that each processor's variables are on their own page */
   //long pad_begin[PAD_SIZE];

   real tnow;         /* current value of simulation time */
   real tout;           /* time next output is due */
   long nstep;        /* number of integration steps so far */

   long workMin, workMax;/* interval of cost to be treated by a proc */

   vector min = {0}, max = {0};   /* min and max of coordinates for each Proc. */

   long mynumcell;  /* num. of cells used for this proc in ctab */
   long mynumleaf;  /* num. of leaves used for this proc in ctab */
   size_t mynbody;    /* num bodies allocated to the processor */
   bodyptr* mybodytab;  /* array of bodies allocated / processor */
   long myncell;  /* num cells allocated to the processor */
   cellptr* mycelltab;  /* array of cellptrs allocated to the processor */
   long mynleaf;  /* number of leaves allocated to the processor */
   leafptr* myleaftab;  /* array of leafptrs allocated to the processor */

   //decltype(celltab)::local_type ctab;  /* array of cells used for the tree. */
   //decltype(leaftab)::local_type ltab;  /* array of cells used for the tree. */

   long myn2bcalc;  /* body-body force calculations for each processor */
   long mynbccalc;  /* body-cell force calculations for each processor */
   long myselfint;  /* count self-interactions for each processor */
   long myn2bterm;  /* count body-body terms for a body */
   long mynbcterm;  /* count body-cell terms for a body */
   bool skipself;   /* true if self-interaction skipped OK */
   bodyptr pskip;       /* body to skip in force evaluation */
   vector pos0 = {0};         /* point at which to evaluate field */
   real phi0 = {0};           /* computed potential at pos0 */
   vector acc0 = {0};         /* computed acceleration at pos0 */
   vector dr = {0};     /* data to be shared */
   real drsq;       /* between gravsub and subdivp */
   nodeptr pmem;  /* remember particle data */

   nodeptr Current_Root; //Cell pointer
   long Root_Coords[NDIM];

   real mymtot;       /* total mass of N-body system */
   real myetot[3];    /* binding, kinetic, potential energy */
   matrix myketen;    /* kinetic energy tensor */
   matrix mypeten;    /* potential energy tensor */
   vector mycmphase[2]; /* center of mass coordinates */
   vector myamvec = {0};    /* angular momentum vector */

   long pad_end[PAD_SIZE];
};

extern struct local_memory Local;

void SlaveStart(void);
void stepsystem(long ProcessId);
void ComputeForces(long ProcessId);
void Help(void);
void ANLinit(void);
void init_root(void);
void tab_init(void);
void startrun(void);
void testdata(void);
void pickshell(real * vec, real rad);
void find_my_initial_bodies(dash::Array<body> & btab, long nbody, long ProcessId);
void find_my_bodies(nodeptr mycell, long work, long direction, long ProcessId);
void Housekeep(long ProcessId);
void setbound(void);
long  Log_base_2(long number);
long  intpow(long i, long j);
#define LOG_MESSAGE(...) DASH_LOG_DEBUG("DASH-BARNES: " __VA_ARGS__);



#endif
