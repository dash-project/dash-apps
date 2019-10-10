/*--------------------------------------------------------------------

    Information on NAS Parallel Benchmarks is available at:

    http://www.nas.nasa.gov/Software/NPB/

    Authors: M. Yarrow
         H. Jin

    CPP and OpenMP version:
            Dalvan Griebler <dalvangriebler@gmail.com>
            Júnior Löff <loffjh@gmail.com>

--------------------------------------------------------------------*/

#include <libdash.h>

#include "npbparams.hpp"
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <thread>


/*****************************************************************/
/* For serial IS, buckets are not really req'd to solve NPB1 IS  */
/* spec, but their use on some machines improves performance, on */
/* other machines the use of buckets compromises performance,    */
/* probably because it is extra computation which is not req'd.  */
/* (Note: Mechanism not understood, probably cache related)      */
/* Example:  SP2-66MhzWN:  50% speedup with buckets              */
/* Example:  SGI Indy5000: 50% slowdown with buckets             */
/* Example:  SGI O2000:   400% slowdown with buckets (Wow!)      */
/*****************************************************************/
/* To disable the use of buckets, comment out the following line */
//#define USE_BUCKETS

/* Uncomment below for cyclic schedule */
/*#define SCHED_CYCLIC*/


/******************/
/* default values */
/******************/
#ifndef CLASS
#define CLASS 'S'
#endif


/*************/
/*  CLASS S  */
/*************/
#if CLASS == 'S'
#define  TOTAL_KEYS_LOG_2    16
#define  MAX_KEY_LOG_2       11
#define  NUM_BUCKETS_LOG_2   9
#endif


/*************/
/*  CLASS W  */
/*************/
#if CLASS == 'W'
#define  TOTAL_KEYS_LOG_2    20
#define  MAX_KEY_LOG_2       16
#define  NUM_BUCKETS_LOG_2   10
#endif

/*************/
/*  CLASS A  */
/*************/
#if CLASS == 'A'
#define  TOTAL_KEYS_LOG_2    23
#define  MAX_KEY_LOG_2       19
#define  NUM_BUCKETS_LOG_2   10
#endif


/*************/
/*  CLASS B  */
/*************/
#if CLASS == 'B'
#define  TOTAL_KEYS_LOG_2    25
#define  MAX_KEY_LOG_2       21
#define  NUM_BUCKETS_LOG_2   10
#endif


/*************/
/*  CLASS C  */
/*************/
#if CLASS == 'C'
#define  TOTAL_KEYS_LOG_2    27
#define  MAX_KEY_LOG_2       23
#define  NUM_BUCKETS_LOG_2   10
#endif


/*************/
/*  CLASS D  */
/*************/
#if CLASS == 'D'
#define  TOTAL_KEYS_LOG_2    31
#define  MAX_KEY_LOG_2       27
#define  NUM_BUCKETS_LOG_2   10
#endif


#if CLASS == 'D'
#define  TOTAL_KEYS          (1L << TOTAL_KEYS_LOG_2)
#else
#define  TOTAL_KEYS          (1 << TOTAL_KEYS_LOG_2)
#endif
#define  MAX_KEY             (1 << MAX_KEY_LOG_2)
#define  NUM_BUCKETS         (1 << NUM_BUCKETS_LOG_2)
#define  NUM_KEYS            TOTAL_KEYS
#define  SIZE_OF_BUFFERS     NUM_KEYS


#define  MAX_ITERATIONS      10
#define  TEST_ARRAY_SIZE     5


/*************************************/
/* Typedef: if necessary, change the */
/* size of int here by changing the  */
/* int type to, say, long            */
/*************************************/
#if CLASS == 'D'
typedef  long INT_TYPE;
#else
typedef  int  INT_TYPE;
#endif


/********************/
/* Some global info */
/********************/
int passed_verification;


/************************************/
/* These are the three main arrays. */
/* See SIZE_OF_BUFFERS def above    */
/************************************/
dash::Array<INT_TYPE> key_array;
dash::Array<INT_TYPE> key_buff2;

INT_TYPE partial_verify_vals[TEST_ARRAY_SIZE];

#ifdef USE_BUCKETS
dash::Array<INT_TYPE> bucket_size;

using pattern_t = dash::CSRPattern<1>;
using array_t   = dash::Array<INT_TYPE, pattern_t::index_type, pattern_t>;

array_t key_buff1;
#else
dash::Array<INT_TYPE> key_buff1;

dash::Array<INT_TYPE> work_buffers;
#endif


/**********************/
/* Partial verif info */
/**********************/
INT_TYPE test_index_array[TEST_ARRAY_SIZE],
         test_rank_array[TEST_ARRAY_SIZE],

         S_test_index_array[TEST_ARRAY_SIZE] =
{48427,17148,23627,62548,4431},
S_test_rank_array[TEST_ARRAY_SIZE] =
{0,18,346,64917,65463},

W_test_index_array[TEST_ARRAY_SIZE] =
{357773,934767,875723,898999,404505},
W_test_rank_array[TEST_ARRAY_SIZE] =
{1249,11698,1039987,1043896,1048018},

A_test_index_array[TEST_ARRAY_SIZE] =
{2112377,662041,5336171,3642833,4250760},
A_test_rank_array[TEST_ARRAY_SIZE] =
{104,17523,123928,8288932,8388264},

B_test_index_array[TEST_ARRAY_SIZE] =
{41869,812306,5102857,18232239,26860214},
B_test_rank_array[TEST_ARRAY_SIZE] =
{33422937,10244,59149,33135281,99},

C_test_index_array[TEST_ARRAY_SIZE] =
{44172927,72999161,74326391,129606274,21736814},
C_test_rank_array[TEST_ARRAY_SIZE] =
{61147,882988,266290,133997595,133525895},

D_test_index_array[TEST_ARRAY_SIZE] =
{1317351170,995930646,1157283250,1503301535,1453734525},
D_test_rank_array[TEST_ARRAY_SIZE] =
{1,36538729,1978098519,2145192618,2147425337};


/***********************/
/* function prototypes */
/***********************/
double	randlc( double *X, double *A );

void full_verify( void );

/*void c_print_results( char   *name,
                      char   class,
                      int    n1,
                      int    n2,
                      int    n3,
                      int    niter,
                      double t,
                      double mops,
		      char   *optype,
                      int    passed_verification,
                      char   *npbversion,
                      char   *compiletime,
                      char   *cc,
                      char   *clink,
                      char   *c_lib,
                      char   *c_inc,
                      char   *cflags,
                      char   *clinkflags );*/
void c_print_results( char   *name, char   class_npb, int    n1, int n2, int n3, int niter, int  nthreads, double t,
                      double mops, char   *optype, int    passed_verification, char   *npbversion, char   *compiletime, char   *cc,
                      char   *clink, char   *c_lib, char   *c_inc, char   *cflags, char   *clinkflags, char   *rand);

void    timer_clear( int n );
void    timer_start( int n );
void    timer_stop( int n );
double  timer_read( int n );


/*
 *    FUNCTION RANDLC (X, A)
 *
 *  This routine returns a uniform pseudorandom double precision number in the
 *  range (0, 1) by using the linear congruential generator
 *
 *  x_{k+1} = a x_k  (mod 2^46)
 *
 *  where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
 *  before repeating.  The argument A is the same as 'a' in the above formula,
 *  and X is the same as x_0.  A and X must be odd double precision integers
 *  in the range (1, 2^46).  The returned value RANDLC is normalized to be
 *  between 0 and 1, i.e. RANDLC = 2^(-46) * x_1.  X is updated to contain
 *  the new seed x_1, so that subsequent calls to RANDLC using the same
 *  arguments will generate a continuous sequence.
 *
 *  This routine should produce the same results on any computer with at least
 *  48 mantissa bits in double precision floating point data.  On Cray systems,
 *  double precision should be disabled.
 *
 *  David H. Bailey     October 26, 1990
 *
 *     IMPLICIT DOUBLE PRECISION (A-H, O-Z)
 *     SAVE KS, R23, R46, T23, T46
 *     DATA KS/0/
 *
 *  If this is the first call to RANDLC, compute R23 = 2 ^ -23, R46 = 2 ^ -46,
 *  T23 = 2 ^ 23, and T46 = 2 ^ 46.  These are computed in loops, rather than
 *  by merely using the ** operator, in order to insure that the results are
 *  exact on all systems.  This code assumes that 0.5D0 is represented exactly.
 */

/*****************************************************************/
/*************           R  A  N  D  L  C             ************/
/*************                                        ************/
/*************    portable random number generator    ************/
/*****************************************************************/

static int      KS=0;
static double	R23, R46, T23, T46;

double	randlc( double *X, double *A )
{
    double		T1, T2, T3, T4;
    double		A1;
    double		A2;
    double		X1;
    double		X2;
    double		Z;
    int     		i, j;

    if (KS == 0)
    {
        R23 = 1.0;
        R46 = 1.0;
        T23 = 1.0;
        T46 = 1.0;

        for (i=1; i<=23; i++)
        {
            R23 = 0.50 * R23;
            T23 = 2.0 * T23;
        }
        for (i=1; i<=46; i++)
        {
            R46 = 0.50 * R46;
            T46 = 2.0 * T46;
        }
        KS = 1;
    }

    /*  Break A into two parts such that A = 2^23 * A1 + A2 and set X = N.  */

    T1 = R23 * *A;
    j  = T1;
    A1 = j;
    A2 = *A - T23 * A1;

    /*  Break X into two parts such that X = 2^23 * X1 + X2, compute
    Z = A1 * X2 + A2 * X1  (mod 2^23), and then
    X = 2^23 * Z + A2 * X2  (mod 2^46).                            */

    T1 = R23 * *X;
    j  = T1;
    X1 = j;
    X2 = *X - T23 * X1;
    T1 = A1 * X2 + A2 * X1;

    j  = R23 * T1;
    T2 = j;
    Z = T1 - T23 * T2;
    T3 = T23 * Z + A2 * X2;
    j  = R46 * T3;
    T4 = j;
    *X = T3 - T46 * T4;
    return(R46 * *X);
}




/*****************************************************************/
/************   F  I  N  D  _  M  Y  _  S  E  E  D    ************/
/************                                         ************/
/************ returns parallel random number seq seed ************/
/*****************************************************************/

/*
 * Create a random number sequence of total length nn residing
 * on np number of processors.  Each processor will therefore have a
 * subsequence of length nn/np.  This routine returns that random
 * number which is the first random number for the subsequence belonging
 * to processor rank kn, and which is used as seed for proc kn ran # gen.
 */

double   find_my_seed( int kn,        /* my processor rank, 0<=kn<=num procs */
                       int np,        /* np = num procs                      */
                       long nn,       /* total num of ran numbers, all procs */
                       double s,      /* Ran num seed, for ex.: 314159265.00 */
                       double a )     /* Ran num gen mult, try 1220703125.00 */
{

    double t1,t2;
    long   mq,nq,kk,ik;

    if ( kn == 0 ) return s;

    mq = (nn/4 + np - 1) / np;
    nq = mq * 4 * kn;               /* number of rans to be skipped */

    t1 = s;
    t2 = a;
    kk = nq;
    while ( kk > 1 ) {
        ik = kk / 2;
        if( 2 * ik ==  kk ) {
            (void)randlc( &t2, &t2 );
            kk = ik;
        }
        else {
            (void)randlc( &t1, &t2 );
            kk = kk - 1;
        }
    }
    (void)randlc( &t1, &t2 );

    return( t1 );

}



/*****************************************************************/
/*************      C  R  E  A  T  E  _  S  E  Q      ************/
/*****************************************************************/

void	create_seq( double seed, double a )
{
    double x, s;
    INT_TYPE i, k;

        INT_TYPE k1, k2;
        double an = a;
        int myid, num_procs;
        INT_TYPE mq;

        myid = dash::myid();
        num_procs = dash::size();

        mq = (NUM_KEYS + num_procs - 1) / num_procs;
        k1 = mq * myid;
        k2 = k1 + mq;
        if ( k2 > NUM_KEYS ) k2 = NUM_KEYS;

        KS = 0;
        s = find_my_seed( myid, num_procs, (long)4*NUM_KEYS, seed, an );

        k = MAX_KEY/4;

        for (i=k1; i<k2; i++) {
            x = randlc(&s, &an);
            x += randlc(&s, &an);
            x += randlc(&s, &an);
            x += randlc(&s, &an);

            key_array[i] = k*x;
        }
}



/*****************************************************************/
/*****************    Allocate Working Buffer     ****************/
/*****************************************************************/
void *alloc_mem( size_t size )
{
    void *p;

    p = (void *)malloc(size);
    if (!p) {
        perror("Memory allocation error");
        exit(1);
    }
    return p;
}

void alloc_key_buff( void )
{
    INT_TYPE i;
    int      num_procs;


    num_procs = dash::size();

		key_array.allocate(NUM_KEYS, dash::BLOCKED);
    key_buff2.allocate(NUM_KEYS, dash::BLOCKED);

#ifdef USE_BUCKETS

		bucket_size.allocate(NUM_BUCKETS*num_procs, dash::BLOCKED);

    int shift = MAX_KEY_LOG_2 - NUM_BUCKETS_LOG_2;
    INT_TYPE num_bucket_keys = (1L << shift);

    std::vector<pattern_t::size_type> l_sizes;
    dash::Array<int> dummy(NUM_BUCKETS);

    dash::Array<pattern_t::size_type> my_elem_count(dash::size());

    my_elem_count.local[0] = dummy.lsize() * num_bucket_keys;

    dash::barrier();

    for (int unit_idx = 0; unit_idx < dash::size(); ++unit_idx) {
      l_sizes.push_back(my_elem_count[unit_idx]);
    }

    pattern_t pattern(l_sizes);

    key_buff1.allocate(pattern);

#else
  key_buff1.allocate(MAX_KEY, dash::BLOCKED);

  work_buffers.allocate(num_procs*MAX_KEY);
#endif /*USE_BUCKETS*/
}



/*****************************************************************/
/*************    F  U  L  L  _  V  E  R  I  F  Y     ************/
/*****************************************************************/


void full_verify( void )
{
    /*  Now, finally, sort the keys:  */

    /*  Copy keys into work array; keys in key_array will be reassigned. */

    std::copy(key_array.lbegin(), key_array.lend(), key_buff2.lbegin());

    key_array.barrier();

    /* This is actual sorting. Each thread is responsible for
    a subset of key values */
    INT_TYPE j = dash::size();
    j = (MAX_KEY + j - 1) / j;
    INT_TYPE k1 = j * dash::myid();
    INT_TYPE k2 = k1 + j;
    if (k2 > MAX_KEY) k2 = MAX_KEY;

    std::for_each(key_buff2.begin(), key_buff2.end(), [k1, k2](INT_TYPE i) {
            if (i >= k1 && i < k2) {
                key_array[--key_buff1[i]] = i;
            }
        });

    dash::barrier();

    /*  Confirm keys correctly sorted: count incorrectly sorted keys, if any */

    dash::Array<INT_TYPE> faults(dash::size());
    dash::fill(faults.begin(), faults.end(), 0);

    dash::Array<INT_TYPE> v2(NUM_KEYS-1);

    //std::iota(v2.lbegin(), v2.lend(), dash::myid() * ceil((double) NUM_KEYS / dash::size())); //only works when using blocking pattern
    if( 0 == dash::myid()) std::iota(v2.begin(), v2.end(), 1);

		v2.barrier();

		dash::for_each(v2.begin(), v2.end(), [&faults, &v2](INT_TYPE i)
		{
        if( key_array[i-1] > key_array[i] ){
          printf("Fault at %d of %d\n", i, (int) v2[NUM_KEYS-1]);
            faults.local[0]++;
          }
    });

    dash::barrier();


    if(0 == dash::myid()) {
      for(int i = 1; i < dash::size(); i++) faults[0] += faults[i];

      if( faults[0] != 0 )
          printf( "Full_verify: number of keys out of sort: %ld\n", (long) faults[0] );
          else
          passed_verification++;
    }
}




/*****************************************************************/
/*************             R  A  N  K             ****************/
/*****************************************************************/


void rank( int iteration )
{

  INT_TYPE    i, k;

#ifdef USE_BUCKETS
    int shift = MAX_KEY_LOG_2 - NUM_BUCKETS_LOG_2;
    INT_TYPE num_bucket_keys = (1L << shift);
#endif

    key_array[iteration] = iteration;
    key_array[iteration+MAX_ITERATIONS] = MAX_KEY - iteration;


    /*  Determine where the partial verify test keys are, load into  */
    /*  top of array bucket_size                                     */
    for( i=0; i<TEST_ARRAY_SIZE; i++ )
      partial_verify_vals[i] = key_array[test_index_array[i]];


    int myid = dash::myid();
    int num_procs = dash::size();

    /*  Bucket sort is known to improve cache performance on some   */
    /*  cache based systems.  But the actual performance may depend */
    /*  on cache size, problem size. */
#ifdef USE_BUCKETS

    //INT_TYPE bucket_ptrs_local[NUM_BUCKETS];
    dash::Array<INT_TYPE> bucket_ptrs(NUM_BUCKETS*num_procs);

    /*  Initialize */
    for( i=0; i<NUM_BUCKETS; i++ )
      bucket_size.local[i] = 0;

    /*  Determine the number of keys in each bucket */
  	dash::for_each(key_array.begin(), key_array.end(), [shift](int k) {
  		bucket_size.local[k >> shift]++;
  	});

    /*  Accumulative bucket sizes are the bucket pointers.
    These are global sizes accumulated upon to each bucket */
    bucket_ptrs.local[0] = 0;
    for( k=0; k< myid; k++ )  {
      bucket_ptrs.local[0] += bucket_size[k*NUM_BUCKETS];
    }

    for( i=1; i< NUM_BUCKETS; i++ ) {
      bucket_ptrs.local[i] = bucket_ptrs.local[i-1];
      for( k=0; k< myid; k++ )
        bucket_ptrs.local[i] += bucket_size[k*NUM_BUCKETS+i];
        for( k=myid; k< num_procs; k++ )
          bucket_ptrs.local[i] += bucket_size[k*NUM_BUCKETS+i-1];
    }
    ////////////////////////////////////
    dash::barrier();
    typedef dash::CSRPattern<1>           pattern_t;
    typedef int                           index_t;
  	typedef typename pattern_t::size_type extent_t;

  	std::vector<extent_t> local_sizes;
    dash::Array<int> my_elem_count(dash::size());
    int local_buckets = ceil((double) NUM_BUCKETS / dash::size());

    my_elem_count.local[0] = ((myid == num_procs - 1)? NUM_KEYS : bucket_ptrs[(myid + 1)*local_buckets]) - bucket_ptrs[(myid + 0)*local_buckets];

    dash::barrier();

    for (int unit_idx = 0; unit_idx < dash::size(); ++unit_idx) {
      local_sizes.push_back(my_elem_count[unit_idx]);
    }

  	pattern_t pattern(local_sizes);

  	dash::Array<INT_TYPE, index_t, pattern_t> key_buff2l(pattern);

    ////////////////////////////////7
  	dash::for_each(key_array.begin(), key_array.end(), [&shift, &bucket_ptrs, &key_buff2l](int k) {
  		key_buff2l[bucket_ptrs.local[k >> shift]++] = k;
  	});

  	if (myid < num_procs-1) {
  		for( i=0; i< NUM_BUCKETS; i++ )
  			for( k=myid+1; k< num_procs; k++ )
  				bucket_ptrs.local[i] += bucket_size[k*NUM_BUCKETS+i];
  	}


    /*  Now, buckets are sorted.  We only need to sort keys inside
    each bucket, which can be done in parallel.  Because the distribution
    of the number of keys in the buckets is Gaussian, the use of
    a dynamic schedule should improve load balance, thus, performance     */

  	dash::Array<int> v(NUM_BUCKETS);

  	//if( 0 == dash::myid()) std::iota(v.begin(), v.end(), 0);

  	std::iota(v.lbegin(), v.lend(), dash::myid() * local_buckets); //only works when using blocking pattern

  	v.barrier();

  	dash::for_each(v.begin(), v.end(), [&num_bucket_keys, &bucket_ptrs, &key_buff2l, &local_buckets, &myid](int i)	{

      int key_buff1_offset = local_buckets * myid * num_bucket_keys;
  		/*  Clear the work array section associated with each bucket */
  		INT_TYPE j;
  		INT_TYPE k1 = i * num_bucket_keys - key_buff1_offset;
  		INT_TYPE k2 = k1 + num_bucket_keys;
  		for ( j = k1; j < k2; j++ )
  			key_buff1.local[j] = 0;

  		/*  Ranking of all keys occurs in this section:                 */
      int key_buff2_offset = (myid > 0)? bucket_ptrs.local[myid*local_buckets-1] : 0;

  		/*  In this section, the keys themselves are used as their
  		own indexes to determine how many of each there are: their
  		individual population                                       */
  		INT_TYPE m = (i > 0)? bucket_ptrs.local[i-1] : 0;
  		for ( j = m - key_buff2_offset; j < bucket_ptrs.local[i] - key_buff2_offset; j++ )
  			key_buff1.local[key_buff2l.local[j] - key_buff1_offset]++;  /* Now they have individual key population */

  			/*  To obtain ranks of each key, successively add the individual key
  			population, not forgetting to add m, the total of lesser keys,
  			to the first key population                                          */
  			key_buff1.local[k1] += m;
  			for ( j = k1+1; j < k2; j++ )
  				key_buff1.local[j] += key_buff1.local[j-1];
  		});

#else /*USE_BUCKETS*/

    std::fill(work_buffers.lbegin(), work_buffers.lend(), 0);

    // compute the histogram for the local keys
    for(int i=0; i<key_array.lsize(); i++) {
        work_buffers.local[ key_array.local[i] ]++;
      }

    // turn it into a cumulative histogram
    for(int i=0; i<MAX_KEY-1; i++ ) {
      work_buffers.local[i+1] += work_buffers.local[i];
    }

    // compute the offset of this unit's local part in
    // the global key_buff1 array
    auto& pat = key_buff1.pattern();
    int goffs = pat.global(0);

    for(int i=0; i<key_buff1.lsize(); i++ ) {
      key_buff1.local[i] = work_buffers.local[goffs+i];
    }

    dash::barrier();

    for(int unit=1; unit<num_procs; unit++ ) {
      for(int i=0; i<key_buff1.lsize(); i++ ) {
        key_buff1.local[i] += work_buffers[((myid+unit)%num_procs)*MAX_KEY+goffs+i];
      }
    }

    dash::barrier();

#endif /*USE_BUCKETS*/


    /* This is the partial verify test section */
    /* Observe that test_rank_array vals are   */
    /* shifted differently for different cases */

		if(0 == myid) {

	    for( i=0; i<TEST_ARRAY_SIZE; i++ )
	    {
	        k = partial_verify_vals[i];          // test vals were put here
	        if( 0 < k  &&  k <= NUM_KEYS-1 )
	        {
	            INT_TYPE key_rank = key_buff1[k-1];
	            int failed = 0;

	            switch( CLASS )
	            {
	            case 'S':
	                if( i <= 2 ) {
	                    if( key_rank != test_rank_array[i]+iteration )
	                        failed = 1;
	                    else
	                        passed_verification++;
	                } else {
	                    if( key_rank != test_rank_array[i]-iteration )
	                        failed = 1;
	                    else
	                        passed_verification++;
	                }
	                break;
	            case 'W':
	                if( i < 2 ) {
	                    if( key_rank != test_rank_array[i]+(iteration-2) )
	                        failed = 1;
	                    else
	                        passed_verification++;
	                } else {
	                    if( key_rank != test_rank_array[i]-iteration )
	                        failed = 1;
	                    else
	                        passed_verification++;
	                }
	                break;
	            case 'A':
	                if( i <= 2 ) {
	                    if( key_rank != test_rank_array[i]+(iteration-1) )
	                        failed = 1;
	                    else
	                        passed_verification++;
	                } else {
	                    if( key_rank != test_rank_array[i]-(iteration-1) )
	                        failed = 1;
	                    else
	                        passed_verification++;
	                }
	                break;
	            case 'B':
	                if( i == 1 || i == 2 || i == 4 ) {
	                    if( key_rank != test_rank_array[i]+iteration )
	                        failed = 1;
	                    else
	                        passed_verification++;
	                } else {
	                    if( key_rank != test_rank_array[i]-iteration )
	                        failed = 1;
	                    else
	                        passed_verification++;
	                }
	                break;
	            case 'C':
	                if( i <= 2 ) {
	                    if( key_rank != test_rank_array[i]+iteration )
	                        failed = 1;
	                    else
	                        passed_verification++;
	                } else {
	                    if( key_rank != test_rank_array[i]-iteration )
	                        failed = 1;
	                    else
	                        passed_verification++;
	                }
	                break;
	            case 'D':
	                if( i < 2 ) {
	                    if( key_rank != test_rank_array[i]+iteration )
	                        failed = 1;
	                    else
	                        passed_verification++;
	                } else {
	                    if( key_rank != test_rank_array[i]-iteration )
	                        failed = 1;
	                    else
	                        passed_verification++;
	                }
	                break;
	            }
	            if( failed == 1 )
	                printf( "Failed partial verification: "
	                        "iteration %d, test key %d\n",
	                        iteration, (int)i );
	        }
	    }
	}

}


/*****************************************************************/
/*************             M  A  I  N             ****************/
/*****************************************************************/

int main( int argc, char **argv )
{
		dash::init(&argc, &argv);

    int nthreads=dash::size();
    int   i, iteration, timer_on;
    double  timecounter;

		if ( 0 == dash::myid() ) {

    	FILE *fp;


    	/*  Initialize timers  */
    	timer_on = 0;
    	if ((fp = fopen("timer.flag", "r")) != NULL) {
        	fclose(fp);
        	timer_on = 1;
    	}
	    timer_clear( 0 );
	    if (timer_on) {
	        timer_clear( 1 );
	        timer_clear( 2 );
	        timer_clear( 3 );
	    }

    	if (timer_on) timer_start( 3 );


	    /*  Initialize the verification arrays if a valid class */
	    for( i=0; i<TEST_ARRAY_SIZE; i++ )
	        switch( CLASS )
	        {
	        case 'S':
	            test_index_array[i] = S_test_index_array[i];
	            test_rank_array[i]  = S_test_rank_array[i];
	            break;
	        case 'A':
	            test_index_array[i] = A_test_index_array[i];
	            test_rank_array[i]  = A_test_rank_array[i];
	            break;
	        case 'W':
	            test_index_array[i] = W_test_index_array[i];
	            test_rank_array[i]  = W_test_rank_array[i];
	            break;
	        case 'B':
	            test_index_array[i] = B_test_index_array[i];
	            test_rank_array[i]  = B_test_rank_array[i];
	            break;
	        case 'C':
	            test_index_array[i] = C_test_index_array[i];
	            test_rank_array[i]  = C_test_rank_array[i];
	            break;
	        case 'D':
	            test_index_array[i] = D_test_index_array[i];
	            test_rank_array[i]  = D_test_rank_array[i];
	            break;
	        };



	    /*  Printout initial NPB info */
	    printf  ( "\n\n NAS Parallel Benchmarks 4.0 OpenMP C++ version - IS Benchmark\n\n" );
	    printf("\n\n Developed by: Dalvan Griebler <dalvan.griebler@acad.pucrs.br>\n");
	    printf( " Size:  %ld  (class %c)\n", (long)TOTAL_KEYS, CLASS );
	    printf( " Iterations:  %d\n", MAX_ITERATIONS );
	    printf( " Number of available threads:  %d\n", std::thread::hardware_concurrency() );
	    printf( "\n" );

			#ifdef USE_BUCKETS
				printf(" This Version of the Benchmark DOES USE Buckets\n\n");
			#else
				printf(" This Version of the Benchmark DOES NOT USE Buckets\n\n");
			#endif

	    if (timer_on) timer_start( 1 );
		}

		alloc_key_buff();

    /*  Generate random number sequence and subsequent keys on all procs */
    create_seq( 314159265.00,                    /* Random number gen seed */
                1220703125.00 );                 /* Random number gen mult */

		if ( 0 == dash::myid() ) {
    	if (timer_on) timer_stop( 1 );
		}


    /*  Do one interation for free (i.e., untimed) to guarantee initialization of
    all data and code pages and respective tables */
    //rank( 1 );

    /*  Start verification counter */
    passed_verification = 0;

    if (0 == dash::myid()) if( CLASS != 'S' ) printf( "\n   iteration\n" );

    /*  Start timer  */
    if (0 == dash::myid()) timer_start( 0 );

    /*  This is the main iteration */
    for( iteration=1; iteration<=MAX_ITERATIONS; iteration++ )
    {
        if (0 == dash::myid()) if( CLASS != 'S' ) printf( "        %d\n", iteration );
        rank( iteration );
        dash::barrier();
    }

		if (0 == dash::myid()) {

	    /*  End of timing, obtain maximum time of all processors */
	    timer_stop( 0 );

	    timecounter = timer_read( 0 );


	    /*  This tests that keys are in sequence: sorting of last ranked key seq
	    occurs here, but is an untimed operation                             */
	    if (timer_on) timer_start( 2 );
		}

	  full_verify();

		if (0 == dash::myid()) {
	    if (timer_on) timer_stop( 2 );

	    if (timer_on) timer_stop( 3 );

	    /*  The final printout  */
	    if( passed_verification != 5*MAX_ITERATIONS + 1 )
	        passed_verification = 0;

	    /*c_print_results( "IS", CLASS, (int)(TOTAL_KEYS/64), 64, 0, MAX_ITERATIONS, timecounter, ((double) (MAX_ITERATIONS*TOTAL_KEYS))
	    /timecounter/1000000., "keys ranked", passed_verification, NPBVERSION, COMPILETIME, CC, CLINK, C_LIB, C_INC,
	    CFLAGS, CLINKFLAGS );*/
	    c_print_results( (char*)"IS", CLASS, TOTAL_KEYS, 0, 0, MAX_ITERATIONS, nthreads, timecounter,
	                     ((double) (MAX_ITERATIONS*TOTAL_KEYS))/timecounter/1000000.0, (char*)"keys ranked", passed_verification,
	                     (char*)NPBVERSION, (char*)COMPILETIME, (char*)CC, (char*)CLINK, (char*)C_LIB, (char*)C_INC, (char*)CFLAGS, (char*)CLINKFLAGS, (char*)"randlc");

	    /*  Print additional timers  */
	    if (timer_on) {
	        double t_total, t_percent;

	        t_total = timer_read( 3 );
	        printf("\nAdditional timers -\n");
	        printf(" Total execution: %8.3f\n", t_total);
	        if (t_total == 0.0) t_total = 1.0;
	        timecounter = timer_read(1);
	        t_percent = timecounter/t_total * 100.;
	        printf(" Initialization : %8.3f (%5.2f%%)\n", timecounter, t_percent);
	        timecounter = timer_read(0);
	        t_percent = timecounter/t_total * 100.;
	        printf(" Benchmarking   : %8.3f (%5.2f%%)\n", timecounter, t_percent);
	        timecounter = timer_read(2);
	        t_percent = timecounter/t_total * 100.;
	        printf(" Sorting        : %8.3f (%5.2f%%)\n", timecounter, t_percent);
	    }

		}

		dash::finalize();

    return 0;
}
/**************************/
/*  E N D  P R O G R A M  */
/**************************/
