#include "npbparams.h"
#include <complex>
#include <string>
#include <libdash.h>

// If processor array is 1x1 -> 0D grid decomposition

// Cache blocking params. These values are good for most
// RISC processors.
// FFT parameters:
//  fftblock controls how many ffts are done at a time.
//  The default is appropriate for most cache-based machines
//  On vector machines, the FFT can be vectorized with vector
//  length equal to the block size, so the block size should
//  be as large as possible. This is the size of the smallest
//  dimension of the problem: 128 for class A, 256 for class B and
//  512 for class C.

#define FFTBLOCK_DEFAULT      32
#define FFTBLOCKPAD_DEFAULT   33

using dcomplex = std::complex<double>;
using std::string;

/* common /blockinfo/ */
static int fftblock, fftblockpad;

// we need a bunch of logic to keep track of how
// arrays are laid out.


// Note: this serial version is the derived from the parallel 0D case
// of the ft NPB.
// The computation proceeds logically as

// set up initial conditions
// fftx(1)
// transpose (1->2)
// ffty(2)
// transpose (2->3)
// fftz(3)
// time evolution
// fftz(3)
// transpose (3->2)
// ffty(2)
// transpose (2->1)
// fftx(1)
// compute residual(1)

// for the 0D, 1D, 2D strategies, the layouts look like xxx
//
//            0D        1D        2D
// 1:        xyz       xyz       xyz

// the array dimensions are stored in dims(coord, phase)
/* common /layout/ */
int dims[3];

#define T_total       1
#define T_setup       2
#define T_fft         3
#define T_evolve      4
#define T_checksum    5
#define T_fftx        6
#define T_ffty        7
#define T_fftz        8
#define T_max         8


// other stuff
/* common /dbg/ */
static bool timers_enabled;
static bool debug;
//static bool debugsynch;

#define SEED          314159265.0
#define A             1220703125.0
#define PI            3.141592653589793238
#define ALPHA         1.0e-6


// roots of unity array
// relies on x being largest dimension?
/* common /ucomm/ */
dash::Array<dcomplex> u;

// for checksum data
/* common /sumcomm/ */
static dcomplex sums[NITER_DEFAULT+1];

// number of iterations
/* common /iter/ */
static int niter;
