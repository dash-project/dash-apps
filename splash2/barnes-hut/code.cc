/*s***********************************************************************/ /*                                                                       */
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
Usage: BARNES <options> < inputfile

Command line options:

    -h : Print out input file description

    Input parameters should be placed in a file and redirected through
    standard input.  There are a total of twelve parameters, and all of
    them have default values.

    1) infile (char*) : The name of an input file that contains particle
       data.

       The format of the file is:
         a) An int representing the number of particles in the distribution
         b) An int representing the dimensionality of the problem (3-D)
         c) A double representing the current time of the simulation
         d) Doubles representing the masses of all the particles
         e) A vector (length equal to the dimensionality) of doubles
            representing the positions of all the particles
         f) A vector (length equal to the dimensionality) of doubles
            representing the velocities of all the particles

       Each of these numbers can be separated by any amount of whitespace.
    2) nbody (int) : If no input file is specified (the first line is
       blank), this number specifies the number of particles to generate
       under a plummer model.  Default is 16384.
    3) seed (int) : The seed used by the random number generator.
       Default is 123.
    4) outfile (char*) : The name of the file that snapshots will be
       printed to. This feature has been disabled in the SPLASH release.
       Default is NULL.
    5) dtime (double) : The integration time-step.
       Default is 0.025.
    6) eps (double) : The usual potential softening
       Default is 0.05.
    7) tol (double) : The cell subdivision tolerance.
       Default is 1.0.
    8) fcells (double) : Number of cells created = fcells * number of
       leaves.
       Default is 2.0.
    9) fleaves (double) : Number of leaves created = fleaves * nbody.
       Default is 0.5.
    10) tstop (double) : The time to stop integration.
       Default is 0.075.
    11) dtout (double) : The data-output interval.
       Default is 0.25.
    12) NPROC (int) : The number of processors.
       Default is 1.
*/

#include <libdash.h>
#include <malloc.h>
#include <pthread.h>
#include <sys/time.h>
#include <unistd.h>
#include <cstdlib>

#include <algorithm>

#include "code.h"
#include "defs.h"
#include "getparam.h"
#include "grav.h"
#include "load.h"
#include "stdinc.h"
#include "util.h"
#include "Logging.h"

typedef struct dash::util::Timer<dash::util::TimeMeasure::Clock> Timer;

/***********************************************
 * DEFINITIONS OF EXTERN GLOBAL VARIABLES      *
 ***********************************************/
string headline;
size_t nbody;
real   dtime;   /* timestep for leapfrog integrator */
real   dtout;   /* time between data outputs */
real   tstop;   /* time to stop calculation */
real   fcells;  /* ratio of cells/leaves allocated */
real   fleaves; /* ratio of leaves/bodies allocated */
real   tol;     /* accuracy parameter: 0.0 => exact */
real   tolsq;   /* square of previous */
real   eps;     /* potential softening parameter */
real   epssq;   /* square of previous */
real   dthf;    /* half time step */

long   maxcell;
long   maxleaf;
size_t maxmybody;
size_t maxmycell;
size_t maxmyleaf;
// dash::Shared<long> n2bcalc; /* total number of body/cell interactions  */
// dash::Shared<long> nbccalc; /* total number of body/body interactions  */
// dash::Shared<long> selfint; /* number of self interactions             */

dash::Array<body> bodytab; /* array size is exactly nbody bodies */
dash::Array<cell> celltab;
dash::Array<leaf> leaftab;

std::vector<dash::Mutex> CellLock;
dash::Mutex              CountLock;

dash::Shared<real> mtot; /* total mass of N-body system             */
// dash::Shared<real> etot[3];      /* binding, kinetic, potential energy */
// dash::Shared<sh_mat> keten; /* kinetic energy tensor                   */
// dash::Shared<sh_mat> peten; /* potential energy tensor                 */
//* center of mass coordinates and velocity */
// dash::Shared<vector> cmphase[2];
// dash::Shared<sh_vec>  amvec;  /* angular momentum vector                 */

dash::Shared<cellptr> G_root; /* root of the whole tree */
/* lower-left corner of coordinate box */
dash::SharedArray<vector> rmin;
/* temporary lower-left corner of the box  */
dash::SharedArray<vector> min;
/* temporary upper right corner of the box */
dash::SharedArray<vector> max;
/* side-length of integer coordinate box   */
dash::Shared<real> rsize;

long computestart, computeend;
long tracktime;
long partitiontime;
long treebuildtime;
long forcecalctime;

struct local_memory Local;

/***********************************************
 * PRIVATE VARIABLES                           *
 ***********************************************/
// static dash::Shared<double> globtout;
// static dash::Shared<double>               globtnow;
// static dash::Shared<int>    globnstep;

const char *defv[] = {
    /* DEFAULT PARAMETER VALUES              */
    /* file names for input/output                                         */
    "in=",  /* snapshot of initial conditions        */
    "out=", /* stream of output snapshots            */

    /* params, used if no input specified, to make a Plummer Model         */
    "nbody=16384", /* number of particles to generate       */
    "seed=123",    /* random number generator seed          */

    /* params to control N-body integration                                */
    "dtime=0.025", /* integration time-step                 */
    "eps=0.05",    /* usual potential softening             */
    "tol=1.0",     /* cell subdivision tolerence            */
    "fcells=2.0",  /* cell allocation parameter             */
    "fleaves=0.5", /* leaf allocation parameter             */

    "tstop=0.075", /* time to stop integration              */
    "dtout=0.25",  /* data-output interval                  */

    "NPROC=1", /* number of processors                  */
};

/* The more complicated 3D case */
#define NUM_DIRECTIONS 32
#define BRC_FUC 0
#define BRC_FRA 1
#define BRA_FDA 2
#define BRA_FRC 3
#define BLC_FDC 4
#define BLC_FLA 5
#define BLA_FUA 6
#define BLA_FLC 7
#define BUC_FUA 8
#define BUC_FLC 9
#define BUA_FUC 10
#define BUA_FRA 11
#define BDC_FDA 12
#define BDC_FRC 13
#define BDA_FDC 14
#define BDA_FLA 15

#define FRC_BUC 16
#define FRC_BRA 17
#define FRA_BDA 18
#define FRA_BRC 19
#define FLC_BDC 20
#define FLC_BLA 21
#define FLA_BUA 22
#define FLA_BLC 23
#define FUC_BUA 24
#define FUC_BLC 25
#define FUA_BUC 26
#define FUA_BRA 27
#define FDC_BDA 28
#define FDC_BRC 29
#define FDA_BDC 30
#define FDA_BLA 31

static long Child_Sequence[NUM_DIRECTIONS][NSUB] = {
    {2, 5, 6, 1, 0, 3, 4, 7}, /* BRC_FUC */
    {2, 5, 6, 1, 0, 7, 4, 3}, /* BRC_FRA */
    {1, 6, 5, 2, 3, 0, 7, 4}, /* BRA_FDA */
    {1, 6, 5, 2, 3, 4, 7, 0}, /* BRA_FRC */
    {6, 1, 2, 5, 4, 7, 0, 3}, /* BLC_FDC */
    {6, 1, 2, 5, 4, 3, 0, 7}, /* BLC_FLA */
    {5, 2, 1, 6, 7, 4, 3, 0}, /* BLA_FUA */
    {5, 2, 1, 6, 7, 0, 3, 4}, /* BLA_FLC */
    {1, 2, 5, 6, 7, 4, 3, 0}, /* BUC_FUA */
    {1, 2, 5, 6, 7, 0, 3, 4}, /* BUC_FLC */
    {6, 5, 2, 1, 0, 3, 4, 7}, /* BUA_FUC */
    {6, 5, 2, 1, 0, 7, 4, 3}, /* BUA_FRA */
    {5, 6, 1, 2, 3, 0, 7, 4}, /* BDC_FDA */
    {5, 6, 1, 2, 3, 4, 7, 0}, /* BDC_FRC */
    {2, 1, 6, 5, 4, 7, 0, 3}, /* BDA_FDC */
    {2, 1, 6, 5, 4, 3, 0, 7}, /* BDA_FLA */

    {3, 4, 7, 0, 1, 2, 5, 6}, /* FRC_BUC */
    {3, 4, 7, 0, 1, 6, 5, 2}, /* FRC_BRA */
    {0, 7, 4, 3, 2, 1, 6, 5}, /* FRA_BDA */
    {0, 7, 4, 3, 2, 5, 6, 1}, /* FRA_BRC */
    {7, 0, 3, 4, 5, 6, 1, 2}, /* FLC_BDC */
    {7, 0, 3, 4, 5, 2, 1, 6}, /* FLC_BLA */
    {4, 3, 0, 7, 6, 5, 2, 1}, /* FLA_BUA */
    {4, 3, 0, 7, 6, 1, 2, 5}, /* FLA_BLC */
    {0, 3, 4, 7, 6, 5, 2, 1}, /* FUC_BUA */
    {0, 3, 4, 7, 6, 1, 2, 5}, /* FUC_BLC */
    {7, 4, 3, 0, 1, 2, 5, 6}, /* FUA_BUC */
    {7, 4, 3, 0, 1, 6, 5, 2}, /* FUA_BRA */
    {4, 7, 0, 3, 2, 1, 6, 5}, /* FDC_BDA */
    {4, 7, 0, 3, 2, 5, 6, 1}, /* FDC_BRC */
    {3, 0, 7, 4, 5, 6, 1, 2}, /* FDA_BDC */
    {3, 0, 7, 4, 5, 2, 1, 6}, /* FDA_BLA */
};

static long Direction_Sequence[NUM_DIRECTIONS][NSUB] = {
    {FRC_BUC, BRA_FRC, FDA_BDC, BLA_FUA, BUC_FLC, FUA_BUC, BRA_FRC, FDA_BLA},
    /* BRC_FUC */
    {FRC_BUC, BRA_FRC, FDA_BDC, BLA_FUA, BRA_FDA, FRC_BRA, BUC_FUA, FLC_BDC},
    /* BRC_FRA */
    {FRA_BDA, BRC_FRA, FUC_BUA, BLC_FDC, BDA_FLA, FDC_BDA, BRC_FRA, FUC_BLC},
    /* BRA_FDA */
    {FRA_BDA, BRC_FRA, FUC_BUA, BLC_FDC, BUC_FLC, FUA_BUC, BRA_FRC, FDA_BLA},
    /* BRA_FRC */
    {FLC_BDC, BLA_FLC, FUA_BUC, BRA_FDA, BDC_FRC, FDA_BDC, BLA_FLC, FUA_BRA},
    /* BLC_FDC */
    {FLC_BDC, BLA_FLC, FUA_BUC, BRA_FDA, BLA_FUA, FLC_BLA, BDC_FDA, FRC_BUC},
    /* BLC_FLA */
    {FLA_BUA, BLC_FLA, FDC_BDA, BRC_FUC, BUA_FRA, FUC_BUA, BLC_FLA, FDC_BRC},
    /* BLA_FUA */
    {FLA_BUA, BLC_FLA, FDC_BDA, BRC_FUC, BLC_FDC, FLA_BLC, BUA_FUC, FRA_BDA},
    /* BLA_FLC */
    {FUC_BLC, BUA_FUC, FRA_BRC, BDA_FLA, BUA_FRA, FUC_BUA, BLC_FLA, FDC_BRC},
    /* BUC_FUA */
    {FUC_BLC, BUA_FUC, FRA_BRC, BDA_FLA, BLC_FDC, FLA_BLC, BUA_FUC, FRA_BDA},
    /* BUC_FLC */
    {FUA_BRA, BUC_FUA, FLC_BLA, BDC_FRC, BUC_FLC, FUA_BUC, BRA_FRC, FDA_BLA},
    /* BUA_FUC */
    {FUA_BRA, BUC_FUA, FLC_BLA, BDC_FRC, BRA_FDA, FRC_BRA, BUC_FUA, FLC_BDC},
    /* BUA_FRA */
    {FDC_BRC, BDA_FDC, FLA_BLC, BUA_FRA, BDA_FLA, FDC_BDA, BRC_FRA, FUC_BLC},
    /* BDC_FDA */
    {FDC_BRC, BDA_FDC, FLA_BLC, BUA_FRA, BUC_FLC, FUA_BUC, BRA_FRC, FDA_BLA},
    /* BDC_FRC */
    {FDA_BLA, BDC_FDA, FRC_BRA, BUC_FLC, BDC_FRC, FDA_BDC, BLA_FLC, FUA_BRA},
    /* BDA_FDC */
    {FDA_BLA, BDC_FDA, FRC_BRA, BUC_FLC, BLA_FUA, FLC_BLA, BDC_FDA, FRC_BUC},
    /* BDA_FLA */

    {BUC_FLC, FUA_BUC, BRA_FRC, FDA_BLA, FUC_BLC, BUA_FUC, FRA_BRC, BDA_FLA},
    /* FRC_BUC */
    {BUC_FLC, FUA_BUC, BRA_FRC, FDA_BLA, FRA_BDA, BRC_FRA, FUC_BUA, BLC_FDC},
    /* FRC_BRA */
    {BRA_FDA, FRC_BRA, BUC_FUA, FLC_BDC, FDA_BLA, BDC_FDA, FRC_BRA, BUC_FLC},
    /* FRA_BDA */
    {BRA_FDA, FRC_BRA, BUC_FUA, FLC_BDC, FRC_BUC, BRA_FRC, FDA_BDC, BLA_FUA},
    /* FRA_BRC */
    {BLC_FDC, FLA_BLC, BUA_FUC, FRA_BDA, FDC_BRC, BDA_FDC, FLA_BLC, BUA_FRA},
    /* FLC_BDC */
    {BLC_FDC, FLA_BLC, BUA_FUC, FRA_BDA, FLA_BUA, BLC_FLA, FDC_BDA, BRC_FUC},
    /* FLC_BLA */
    {BLA_FUA, FLC_BLA, BDC_FDA, FRC_BUC, FUA_BRA, BUC_FUA, FLC_BLA, BDC_FRC},
    /* FLA_BUA */
    {BLA_FUA, FLC_BLA, BDC_FDA, FRC_BUC, FLC_BDC, BLA_FLC, FUA_BUC, BRA_FDA},
    /* FLA_BLC */
    {BUC_FLC, FUA_BUC, BRA_FRC, FDA_BLA, FUA_BRA, BUC_FUA, FLC_BLA, BDC_FRC},
    /* FUC_BUA */
    {BUC_FLC, FUA_BUC, BRA_FRC, FDA_BLA, FLC_BDC, BLA_FLC, FUA_BUC, BRA_FDA},
    /* FUC_BLC */
    {BUA_FRA, FUC_BUA, BLC_FLA, FDC_BRC, FUC_BLC, BUA_FUC, FRA_BRC, BDA_FLA},
    /* FUA_BUC */
    {BUA_FRA, FUC_BUA, BLC_FLA, FDC_BRC, FRA_BDA, BRC_FRA, FUC_BUA, BLC_FDC},
    /* FUA_BRA */
    {BDC_FRC, FDA_BDC, BLA_FLC, FUA_BRA, FDA_BLA, BDC_FDA, FRC_BRA, BUC_FLC},
    /* FDC_BDA */
    {BDC_FRC, FDA_BDC, BLA_FLC, FUA_BRA, FRC_BUC, BRA_FRC, FDA_BDC, BLA_FUA},
    /* FDC_BRC */
    {BDA_FLA, FDC_BDA, BRC_FRA, FUC_BLC, FDC_BRC, BDA_FDC, FLA_BLC, BUA_FRA},
    /* FDA_BDC */
    {BDA_FLA, FDC_BDA, BRC_FRA, FUC_BLC, FLA_BUA, BLC_FLA, FDC_BDA, BRC_FUC},
    /* FDA_BLA */
};

static void LocalInit()
{
  // Command Line Params
  dtime   = 0.025;
  dthf    = 0.5 * dtime;
  dtout   = 0.25;
  tstop   = 0.075;
  fcells  = 2.0;
  fleaves = 5.0;
  tol     = 1.0;
  tolsq   = tol * tol;
  eps     = 0.05;
  epssq   = eps * eps;

  // for allocation of body and leaf tab
  maxleaf = fleaves * nbody;
  maxcell = fcells * maxleaf;

  Local.tnow  = 0.0;
  Local.tout  = Local.tnow + dtout;
  Local.nstep = 0;
}
/***********************************************
 * IMPLEMENTATION OF BARNES HUT                *
 ***********************************************/

int main(int argc, char *argv[])
{
  dash::init(&argc, &argv);

  Timer::Calibrate(0);

  if (argc < 2) {
    if (dash::myid() == 0) {
      std::cout << "usage: " << std::string(argv[0]) << std::endl;
    }
    dash::finalize();
  }

  int  exit = 2;
  long c;

  while ((c = getopt(argc, argv, "h")) != -1) {
    switch (c) {
      case 'h':
        if (dash::myid() == 0) Help();
        exit = EXIT_SUCCESS;
        break;
      default:
        if (dash::myid() == 0) {
          Help();
          std::cerr << "Only valid option is \"-h\".\n";
        }
        exit = EXIT_FAILURE;
        break;
    }
  }

  if (argv[1][0] == '-') {
    if (dash::myid() == 0)
      std::cerr << "argument for nbody is a negative number!\n";

    exit = EXIT_FAILURE;
  }

  if (exit != 2) {
    dash::finalize();
    return exit;
  }
  nbody = ::std::atoi(argv[1]);

  if (argc > 2) {
    int debug = ::std::atoi(argv[2]);

    if (dash::myid() == 0 && debug) {
      int wait = 1;
      while (wait)
        ;
    }
  }

  ANLinit();  // Prepare all global data structures

  // Prepare global data structures
  LocalInit();

  if (dash::myid() == 0) {
    initparam(defv);
    // initoutput();
    startrun();
  }

  tab_init();

  dash::barrier();

  if (dash::myid() == 0) {
    tracktime     = 0;
    partitiontime = 0;
    treebuildtime = 0;
    forcecalctime = 0;
    computestart  = Timer::Now();
    printf(
        "COMPUTESTART  = %12lu\n", static_cast<unsigned long>(computestart));
  }

  SlaveStart();

  // wait for all processes to finish
  dash::barrier();

  if (dash::myid() == 0) {
    computeend = Timer::Now();

    std::ostringstream os;

    os << "COMPUTEEND=" << std::setw(12)
       << static_cast<unsigned long>(computeend) << "\n";
    os << "COMPUTETIME=" << std::setw(12)
       << static_cast<unsigned long>(computeend - computestart) << "\n";

    std::cout << os.str();

    printf("TRACKTIME     = %12lu\n", tracktime);
    printf(
        "PARTITIONTIME = %12lu\t%5.2f\n",
        partitiontime,
        (static_cast<double>(partitiontime) / tracktime));
    printf(
        "TREEBUILDTIME = %12lu\t%5.2f\n",
        treebuildtime,
        (static_cast<double>(treebuildtime) / tracktime));
    printf(
        "FORCECALCTIME = %12lu\t%5.2f\n",
        forcecalctime,
        (static_cast<double>(forcecalctime)) / tracktime);

    printf(
        "RESTTIME      = %12lu\t%5.2f\n",
        tracktime - partitiontime - treebuildtime - forcecalctime,
        static_cast<double>(
            (tracktime - partitiontime - treebuildtime - forcecalctime)) /
            tracktime);
  }
  return EXIT_SUCCESS;
}

/*
 * ANLINIT : initialize ANL macros
 */
void ANLinit()
{
#if 0
  {
    ;
  };
  /* Allocate global, shared memory */

  Global = (struct GlobalMemory *)valloc(sizeof(struct GlobalMemory));
  ;
  if (Global == NULL) error("No initialization for Global\n");

  {
    pthread_barrier_init(&(Global->Barrier), NULL, NPROC);
  };

  {
    pthread_mutex_init(&(Global->CountLock), NULL);
  };
  {
    pthread_mutex_init(&(Global->io_lock), NULL);
  };
#endif

  /*
  for (auto u = 0; u < MAXLOCK; ++u) {
    CellLock.emplace_back(dash::Team::All());
  }
  */
  CellLock.reserve(MAXLOCK);

  std::generate_n(std::back_inserter(CellLock), MAXLOCK, []() {
    return dash::Mutex{dash::Team::All()};
  });

  auto const nprocs = dash::size();

  if (!bodytab.allocate(nbody)) {
    error("testdata: not enough memory\n");
  }

  std::uninitialized_fill(bodytab.lbegin(), bodytab.lend(), body{});

  CountLock.init();

  //(dtime.allocate());
  //(dthf.allocate());
  //(eps.allocate());
  //(epssq.allocate());
  //(tol.allocate());
  //(tolsq.allocate());
  //(fcells.allocate());
  //(fleaves.allocate());
  //(tstop.allocate());
  //(dtout.allocate());
  //(globtout.allocate());
  //(globtnow.allocate());
  //(globnstep.allocate());
  //(maxleaf.allocate());
  //(maxcell.allocate());
  (rsize.init());  //
  (rmin.init());
  (max.init());
  (min.init());
  //(createstart.allocate());
  //(createend.allocate());
  // (computestart.allocate());
  // (computeend.allocate());
  // (trackstart.allocate());
  // (trackend.allocate());
  // (tracktime.allocate());
  // (partitionstart.allocate());
  // (partitionend.allocate());
  // (partitiontime.allocate());
  // (treebuildstart.allocate());
  // (treebuildend.allocate());
  // (treebuildtime.allocate());
  // (forcecalcstart.allocate());
  // (forcecalcend.allocate());
  // (forcecalctime.allocate());
  (G_root.init());
}

/*
 * INIT_ROOT: Processor 0 reinitialize the global root at each time step
 */
void init_root()
{
  // The root is always the first cell
  //
  auto &root  = celltab.local[0];
  root.type   = CELL;
  root.seqnum = 0;
  root.done   = false;
  root.level  = IMAX >> 1;

  //convert reference to pointer
  G_root.set(static_cast<cellptr>(celltab.begin()));
  // The root has initially no children
  std::fill(std::begin(root.subp), std::end(root.subp), cellptr{});
  Local.mynumcell = 1;
}

void tab_init()
{
  ASSERT(maxleaf);
  ASSERT(maxcell);

  celltab.allocate(maxcell);
  leaftab.allocate(maxleaf);

  std::uninitialized_fill(celltab.lbegin(), celltab.lend(), cell{});
  std::uninitialized_fill(leaftab.lbegin(), leaftab.lend(), leaf{});

  /*allocate space for personal lists of body pointers */
  maxmybody = (nbody + maxleaf * MAX_BODIES_PER_LEAF) / bodytab.team().size();

  /* space is allocated so that every */
  /* process can have a maximum of maxmybody pointers to bodies */
  /* then there is an array of bodies called bodytab which is  */
  /* allocated in the distribution generation or when the distr. */
  /* file is read */
  maxmycell = maxcell / bodytab.team().size();
  maxmyleaf = maxleaf / bodytab.team().size();

  Local.mybodytab.resize(maxmybody);
  Local.mycelltab.resize(maxmycell);
  Local.myleaftab.resize(maxmyleaf);

  LOG("maxmybody: " << maxmybody);
  LOG("maxmycell: " << maxmycell);
  LOG("maxmyleaf: " << maxmyleaf);
}

/*
 * SLAVESTART: main task for each processor
 */
void SlaveStart()
{
  long ProcessId = dash::myid();

/* POSSIBLE ENHANCEMENT:  Here is where one might pin processes to
   processors to avoid migration */
#if 0
  /* initialize mybodytabs */
  Local.mybodytab = Local[0].mybodytab + (maxmybody * ProcessId);
  /* note that every process has its own copy   */
  /* of mybodytab, which was initialized to the */
  /* of mybodytab, which was initialized to the */
  /* before create                              */
  Local.mycelltab = Local[0].mycelltab + (maxmycell * ProcessId);
  Local.myleaftab = Local[0].myleaftab + (maxmyleaf * ProcessId);
  /* POSSIBLE ENHANCEMENT:  Here is where one might distribute the
     data across physically distributed memories as desired.

     One way to do this is as follows:

     long i;

     if (ProcessId == 0) {
       for (i=0;i<NPROC;i++) {
         Place all addresses x such that
           &(Local[i]) <= x < &(Local[i])+
             sizeof(struct local_memory) on node i
         Place all addresses x such that
           &(Local[i].mybodytab) <= x < &(Local[i].mybodytab)+
             maxmybody * sizeof(bodyptr) - 1 on node i
         Place all addresses x such that
           &(Local[i].mycelltab) <= x < &(Local[i].mycelltab)+
             maxmycell * sizeof(cellptr) - 1 on node i
         Place all addresses x such that
           &(Local[i].myleaftab) <= x < &(Local[i].myleaftab)+
             maxmyleaf * sizeof(leafptr) - 1 on node i
       }
     }

     barrier(Global->Barstart,NPROC);

  */

#endif
  // Local.tout = globtout.get();
  // Local.tnow  = globtnow.get();
  // Local.nstep = globnstep.get();

  find_my_initial_bodies(bodytab, nbody, ProcessId);

  /* main loop */
  while (Local.tnow < tstop + 0.1 * dtime) {
    stepsystem(ProcessId);
    //if (ProcessId == 0) printtree(G_root.get().get());
  }
  if (ProcessId == 0) printtree(G_root.get().get());
}

void startrun()
{
  long seed;

  // Original Settings as in reference implementation
  seed = 123;
  // outfile        = NULL;

  pranset(seed);
  testdata();
  setbound();

  printf("----------PARAMS----------------\n");
  printf("infile = \n");
  printf("nbody = %ld\n", nbody);
  printf("seed = %ld\n", seed);
  printf("outfile = \n");
  printf("dtime = %f\n", dtime);
  printf("eps = %f\n", eps);
  printf("tol = %f\n", tol);
  printf("fcells = %f\n", fcells);
  printf("fleaves = %f\n", fleaves);
  printf("tstop = %f\n", tstop);
  printf("dtout = %f\n", dtout);
  printf("NPROC = %zd\n", dash::size());
}

/*
 * TESTDATA: generate Plummer model initial conditions for test runs,
 * scaled to units such that M = -4E = G = 1 (Henon, Hegge, etc).
 * See Aarseth, SJ, Henon, M, & Wielen, R (1974) Astr & Ap, 37, 183.
 */

#define MFRAC 0.999 /* mass cut off at MFRAC of total */

void testdata()
{
  real   rsc, vsc, r, v, x, y;
  vector cmr, cmv;
  long   rejects = 0;
  long   halfnbody, i;
  float  offset;

  headline = "Hack code: Plummer model";
  // Local.tnow = 0.0;

  rsc = 9 * PI / 16;
  vsc = sqrt(1.0 / rsc);

  CLRV(cmr);
  CLRV(cmv);

  halfnbody = nbody / 2;
  if (nbody % 2 != 0) halfnbody++;

  for (auto iter = bodytab.begin(); iter < bodytab.begin() + halfnbody;
       ++iter) {
    // local copy
    auto p = static_cast<body>(*iter);

    p.type = BODY;
    // Type = BODY;
    p.mass = 1.0 / nbody;
    p.cost = 1;

    r = 1 / std::sqrt(std::pow(xrand(0.0, MFRAC), -2.0 / 3.0) - 1);
    /*   reject radii greater than 10 */
    while (r > 9.0) {
      rejects++;
      r = 1 / std::sqrt(std::pow(xrand(0.0, MFRAC), -2.0 / 3.0) - 1);
    }
    pickshell(p.pos, rsc * r);
    ADDV(cmr, cmr, p.pos);
    do {
      x = xrand(0.0, 1.0);
      y = xrand(0.0, 0.1);

    } while (y > x * x * std::pow(1 - x * x, 3.5));

    v = sqrt(2.0) * x / std::pow(1 + r * r, 0.25);
    pickshell(p.vel, vsc * v);
    ADDV(cmv, cmv, p.vel);
    // write back
    *iter = p;
  }

  offset = 4.0;

  for (auto iter = bodytab.begin() + halfnbody;
       iter < bodytab.begin() + nbody;
       ++iter) {
    body p = *iter;
    p.type = BODY;
    p.mass = 1.0 / nbody;
    p.cost = 1;
    // p.child_num = 0;

    auto const cp = static_cast<body>(*(iter - halfnbody));
    for (i = 0; i < NDIM; i++) {
      p.pos[i] = cp.pos[i] + offset;
      p.vel[i] = cp.vel[i];
    }
    ADDV(cmr, cmr, p.pos);
    ADDV(cmv, cmv, p.vel);

    *iter = p;
  }

  DIVVS(cmr, cmr, (real)nbody);
  DIVVS(cmv, cmv, (real)nbody);

  for (auto iter = bodytab.begin(); iter < bodytab.end(); ++iter) {
    auto p = static_cast<body>(*iter);
    SUBV(p.pos, p.pos, cmr);
    SUBV(p.vel, p.vel, cmv);
    *iter = p;
  }
}

/*
 * PICKSHELL: pick a random point on a sphere of specified radius.
 */

void pickshell(real *vec, real rad)
{
  long   k;
  double rsq, rsc;

  do {
    for (k = 0; k < NDIM; k++) {
      vec[k] = xrand(-1.0, 1.0);
    }
    DOTVP(rsq, vec, vec);
  } while (rsq > 1.0);

  rsc = rad / sqrt(rsq);
  MULVS(vec, vec, rsc);
}

long intpow(long i, long j)
{
  long k;
  long temp = 1;

  for (k = 0; k < j; k++) temp = temp * i;
  return temp;
}
/*
 * STEPSYSTEM: advance N-body system one time-step.
 */

void stepsystem(long ProcessId)
{
  long i;
  // real Cavg;
  // bodyptr p, *pp;
  vector dvel, vel1, dpos;
  long   trackstart, trackend;
  long   partitionstart, partitionend;
  long   treebuildstart, treebuildend;
  long   forcecalcstart, forcecalcend;

  auto const nprocs = bodytab.team().size();

  if (Local.nstep == 2) {
    /* POSSIBLE ENHANCEMENT:  Here is where one might reset the
       statistics that one is measuring about the parallel execution */
  }

  if ((ProcessId == 0) && (Local.nstep >= 2)) {
    trackstart = Timer::Now();
  }

  if (ProcessId == 0) {
    // Process 0 is initial owner of root
    init_root();
  }
  else {
    Local.mynumcell = 0;
    Local.mynumleaf = 0;
  }

  DASH_LOG_DEBUG("after initializing root in stepsystem");

  /* start at same time */
  /*
  {
    pthread_barrier_wait(&(Global->Barrier));
  };
  */

  dash::barrier();

  /*
  if ((ProcessId == 0) && (Local[ProcessId].nstep >= 2)) {
    {
      struct timeval FullTime;

      gettimeofday(&FullTime, NULL);

      (treebuildstart) =
          (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
    };
  }
  */
  if ((ProcessId == 0) && (Local.nstep >= 2)) {
    treebuildstart = Timer::Now();
  }

  // std::cout << "BEFORE MAKETREE:" << std::endl;
  // printtree(G_root.get().get());
  /* load bodies into tree   */
  maketree(ProcessId);
  // std::cout << "AFTER MAKETREE:" << std::endl;
  // printtree(G_root.get().get());
  // std::cout << cell_ref.get();
  /*
  if ((ProcessId == 0) && (Local[ProcessId].nstep >= 2)) {
    {
      struct timeval FullTime;

      gettimeofday(&FullTime, NULL);

      (treebuildend) =
          (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
    };
    Global->treebuildtime += treebuildend - treebuildstart;
  }
  */
  if (ProcessId == 0 && Local.nstep >= 2) {
    treebuildend = Timer::Now();

    treebuildtime = treebuildend - treebuildstart;
  }

  Housekeep(ProcessId);

  auto const Cavg =
      static_cast<real>(Cost(*(static_cast<cellptr>(G_root.get())))) / nprocs;

  LOG("Cavg: " << Cavg);

  //even distribution
  Local.workMin = Cavg * ProcessId;
  Local.workMax =
      (Cavg * (ProcessId + 1) + (ProcessId == static_cast<long>(nprocs - 1)));

  if ((ProcessId == 0) && (Local.nstep >= 2)) {
    /*
  {
    struct timeval FullTime;

    gettimeofday(&FullTime, NULL);

    (partitionstart) =
        (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
  };
  */
    partitionstart = Timer::Now();
  }

  Local.mynbody = 0;
  find_my_bodies(G_root.get().get(), 0, BRC_FUC, ProcessId);
  LOG("my nbodies: " << Local.mynbody);

  if ((ProcessId == 0) && (Local.nstep >= 2)) {
    /*
    {
      struct timeval FullTime;

      gettimeofday(&FullTime, NULL);

      (partitionend) =
          (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
    };
    */
    partitionend = Timer::Now();

    partitiontime += partitionend - partitionstart;
  }

  if ((ProcessId == 0) && (Local.nstep >= 2)) {
    /*
    {
      struct timeval FullTime;

      gettimeofday(&FullTime, NULL);

      (forcecalcstart) =
          (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
    };
    */
    forcecalcstart = Timer::Now();
  }

  ComputeForces(ProcessId);

  if ((ProcessId == 0) && (Local.nstep >= 2)) {
    /*
    {
      struct timeval FullTime;

      gettimeofday(&FullTime, NULL);

      (forcecalcend) =
          (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
    };
    */
    forcecalcend = Timer::Now();
    forcecalctime += forcecalcend - forcecalcstart;
  }

  /* advance my bodies */
  for (auto pp = std::begin(Local.mybodytab);
       pp < std::next(std::begin(Local.mybodytab), Local.mynbody);
       ++pp) {
    auto p_gptr = *pp;
    auto p_val  = static_cast<body>(*p_gptr);
    MULVS(dvel, p_val.acc, dthf);
    ADDV(vel1, p_val.vel, dvel);
    MULVS(dpos, vel1, dtime);
    ADDV(p_val.pos, p_val.pos, dpos);
    ADDV(p_val.vel, vel1, dvel);

    *p_gptr = p_val;

    for (i = 0; i < NDIM; i++) {
      if (p_val.pos[i] < Local.min[i]) {
        Local.min[i] = p_val.pos[i];
      }
      if (p_val.pos[i] > Local.max[i]) {
        Local.max[i] = p_val.pos[i];
      }
    }
  }

  /*
  {
    pthread_mutex_lock(&(Global->CountLock));
  };
  */

  CountLock.lock();
  // TODO rko: This can be replaced with an atomic operation
  vector g_min_v, g_max_v, g_min_v_cpy, g_max_v_cpy;

  // Get the global minimum into a local variable
  min.get(g_min_v);
  max.get(g_max_v);

  // Keep a copy of it
  std::copy(std::begin(g_min_v), std::end(g_min_v), std::begin(g_min_v_cpy));
  std::copy(std::begin(g_max_v), std::end(g_max_v), std::begin(g_max_v_cpy));

  // Update global min / max
  for (i = 0; i < NDIM; i++) {
    g_min_v[i] = std::min(Local.min[i], g_min_v[i]);
    g_max_v[i] = std::max(Local.max[i], g_max_v[i]);
  }

  // Reset global min and max if value has changed
  if (!(std::equal(
          std::begin(g_min_v), std::end(g_min_v), std::begin(g_min_v_cpy)))) {
    min.set(g_min_v);
  }

  if (!(std::equal(
          std::begin(g_max_v), std::end(g_max_v), std::begin(g_max_v_cpy)))) {
    max.set(g_max_v);
  }



  /*
  {
    pthread_mutex_unlock(&(Global->CountLock));
  };
  */
  CountLock.unlock();

  /* bar needed to make sure that every process has computed its min */
  /* and max coordinates, and has accumulated them into the global   */
  /* min and max, before the new dimensions are computed
   */
  /*
  {
    pthread_barrier_wait(&(Global->Barrier));
  };
  */
  dash::barrier();

  if ((ProcessId == 0) && (Local.nstep >= 2)) {
    /*
    {
      struct timeval FullTime;

      gettimeofday(&FullTime, NULL);

      (trackend) =
          (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
    };
    */
    trackend = Timer::Now();
    tracktime += trackend - trackstart;
  }

  if (ProcessId == 0) {
    // real g_rsize_v = 0;

    // update global max
    max.get(g_max_v);
    min.get(g_min_v);

    // SUBV(Global->max,Global->max,Global->min);
    SUBV(g_max_v, g_max_v, g_min_v);

    /*
    for (i = 0; i < NDIM; i++) {
      if (Global->rsize < Global->max[i]) {
        Global->rsize = Global->max[i];
      }
    }
    */

    auto g_rsize_v = std::max(
        real{0}, *std::max_element(std::begin(g_max_v), std::end(g_max_v)));

    // ADDVS(Global->rmin,Global->min,-Global->rsize/100000.0);
    vector g_rmin_v;
    rmin.get(g_rmin_v);

    // ADDVS(Global->rmin,Global->min,-Global->rsize/100000.0);
    ADDVS(g_rmin_v, g_min_v, -g_rsize_v / 100000.0);
    rmin.set(g_rmin_v);

    // Global->rsize = 1.00002*Global->rsize;
    g_rsize_v = 1.00002 * g_rsize_v;
    rsize.set(g_rsize_v);

    SETVS(g_min_v, 1E99);
    SETVS(g_max_v, -1E99);
    min.set(g_min_v);
    max.set(g_max_v);
  }
  Local.nstep++;
  Local.tnow += dtime;
}

void ComputeForces(long ProcessId)
{
  bodyptr p;
  vector  acc1, dacc, dvel;

  for (auto pp = std::begin(Local.mybodytab);
       pp < std::next(std::begin(Local.mybodytab), Local.mynbody);
       ++pp) {
    p = *pp;
    ASSERT(p != bodyptr{});
    auto p_val = static_cast<body>(*p);
    SETV(acc1, p_val.acc);
    Cost(*p) = 0;
    hackgrav(p, ProcessId);
    // reload p_val since the global pointer p is modified
    p_val = *p;
    Local.myn2bcalc += Local.myn2bterm;
    Local.mynbccalc += Local.mynbcterm;
    if (!Local.skipself) { /*   did we miss self-int?  */
      Local.myselfint++;   /*   count another goofup   */
    }
    if (Local.nstep > 0) {
      /*   use change in accel to make 2nd order correction to vel      */
      SUBV(dacc, p_val.acc, acc1);
      MULVS(dvel, dacc, dthf);
      ADDV(p_val.vel, p_val.vel, dvel);
    }
    *p = p_val;
  }
}

/*
 * FIND_MY_INITIAL_BODIES: puts into mybodytab the initial list of bodies
 * assigned to the processor.
 */

void find_my_initial_bodies(
    dash::Array<body> &btab, long nbody, long ProcessId)
{
  long /*extra,*/ offset;
  size_t          i;

  offset = btab.pattern().global(0);
  // Local.mynbody = nbody / bodytab.team().size();
  // extra         = nbody % bodytab.team().size();
  Local.mynbody = btab.lsize();

#if 0
  if (ProcessId < extra) {
    //Local.mynbody++;
    //offset = Local.mynbody * ProcessId;
  }

  if (ProcessId >= extra) {
    /*
    offset = (Local.mynbody + 1) * extra +
             (ProcessId - extra) * Local.mynbody;
             */
    // local to global index mapping
  }
#endif

  //ASSERT(Local.mynbody <= btab.lsize());

  using iter_t = decltype(btab.begin());


  //TODO: use dash::transform
  for (i = 0; i < Local.mynbody; i++) {
    // copy local pointer
    Local.mybodytab[i] =
        static_cast<typename iter_t::pointer>(btab.begin() + offset + i);
  }
}

void find_my_bodies(nodeptr mycell, long work, long direction, long ProcessId)
{
  long i;

  /* NOTE: Local.mynbody has been reset to 0 before calling this method */
  if (Type((*mycell)) == LEAF) {
    auto       l     = static_cast<leafptr>(mycell);
    auto const l_val = static_cast<leaf>(*l);
    for (i = 0; i < l_val.num_bodies; i++) {
      if (work >= Local.workMin - .1) {
        if ((Local.mynbody + 2) > maxmybody) {
          error(
              "find_my_bodies: Processor %ld needs more than %ld bodies; "
              "increase fleaves\n",
              ProcessId,
              maxmybody);
        }
        Local.mybodytab[Local.mynbody++] = l_val.bodyp[i];
      }
      work += Cost(*(l_val.bodyp[i])).get();
      if (work >= Local.workMax - .1) {
        break;
      }
    }
  }
  else {
    for (i = 0; (i < NSUB) && (work < (Local.workMax - .1)); i++) {
      auto const mycell_val = static_cast<cell>(*NODE_AS_CELL(mycell));

      auto qptr    = mycell_val.subp[Child_Sequence[direction][i]];
      if (qptr != nodeptr{}) {
        auto const cost = static_cast<long>(Cost(*qptr));
        if ((work + cost) >= (Local.workMin - .1)) {
          find_my_bodies(
              qptr, work, Direction_Sequence[direction][i], ProcessId);
        }
        work += static_cast<long>(Cost(*qptr));
      }
    }
  }
}

/*
 * HOUSEKEEP: reinitialize the different variables (in particular global
 * variables) between each time step.
 */

void Housekeep(long ProcessId)
{
  Local.myn2bcalc = Local.mynbccalc = Local.myselfint = 0;
  SETVS(Local.min, 1E99);
  SETVS(Local.max, -1E99);
}

/*
 * SETBOUND: Compute the initial size of the root of the tree; only done
 * before first time step, and only processor 0 does it
 */
void setbound()
{
  long    i;
  real    side;
  bodyptr p;

  SETVS(Local.min, 1E99);
  SETVS(Local.max, -1E99);
  side = 0;

  for (auto iter = bodytab.begin(); iter < bodytab.end(); ++iter) {
    auto const p = static_cast<body>(*iter);
    for (i = 0; i < NDIM; i++) {
      if (p.pos[i] < Local.min[i]) Local.min[i] = p.pos[i];
      if (p.pos[i] > Local.max[i]) Local.max[i] = p.pos[i];
    }
  }

  SUBV(Local.max, Local.max, Local.min);
  for (i = 0; i < NDIM; i++) {
    if (side < Local.max[i]) side = Local.max[i];
  }

  rsize.set(1.00002 * side);

  vector tmp;

  ADDVS(tmp, Local.min, -side / 100000.0);
  rmin.set(tmp);

  SETVS(tmp, -1E99)
  max.set(tmp);

  SETVS(tmp, 1E99)
  min.set(tmp);
}

void Help()
{
  printf(
      "There are a total of twelve parameters, and all of them have default "
      "values.\n");
  printf("\n");
  printf(
      "1) infile (char*) : The name of an input file that contains particle "
      "data.  \n");
  printf("    The format of the file is:\n");
  printf(
      "\ta) An int representing the number of particles in the "
      "distribution\n");
  printf(
      "\tb) An int representing the dimensionality of the problem (3-D)\n");
  printf("\tc) A double representing the current time of the simulation\n");
  printf("\td) Doubles representing the masses of all the particles\n");
  printf("\te) A vector (length equal to the dimensionality) of doubles\n");
  printf("\t   representing the positions of all the particles\n");
  printf("\tf) A vector (length equal to the dimensionality) of doubles\n");
  printf("\t   representing the velocities of all the particles\n");
  printf("\n");
  printf(
      "    Each of these numbers can be separated by any amount of "
      "whitespace.\n");
  printf("\n");
  printf(
      "2) nbody (int) : If no input file is specified (the first line is "
      "blank), this\n");
  printf(
      "    number specifies the number of particles to generate under a "
      "plummer model.\n");
  printf("    Default is 16384.\n");
  printf("\n");
  printf("3) seed (int) : The seed used by the random number generator.\n");
  printf("    Default is 123.\n");
  printf("\n");
  printf(
      "4) outfile (char*) : The name of the file that snapshots will be "
      "printed to. \n");
  printf("    This feature has been disabled in the SPLASH release.\n");
  printf("    Default is NULL.\n");
  printf("\n");
  printf("5) dtime (double) : The integration time-step.\n");
  printf("    Default is 0.025.\n");
  printf("\n");
  printf("6) eps (double) : The usual potential softening\n");
  printf("    Default is 0.05.\n");
  printf("\n");
  printf("7) tol (double) : The cell subdivision tolerance.\n");
  printf("    Default is 1.0.\n");
  printf("\n");
  printf(
      "8) fcells (double) : The total number of cells created is equal to "
      "\n");
  printf("    fcells * number of leaves.\n");
  printf("    Default is 2.0.\n");
  printf("\n");
  printf(
      "9) fleaves (double) : The total number of leaves created is equal to  "
      "\n");
  printf("    fleaves * nbody.\n");
  printf("    Default is 0.5.\n");
  printf("\n");
  printf("10) tstop (double) : The time to stop integration.\n");
  printf("    Default is 0.075.\n");
  printf("\n");
  printf("11) dtout (double) : The data-output interval.\n");
  printf("    Default is 0.25.\n");
  printf("\n");
  printf("12) NPROC (int) : The number of processors.\n");
  printf("    Default is 1.\n");
}
