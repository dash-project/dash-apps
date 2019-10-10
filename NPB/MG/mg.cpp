/*--------------------------------------------------------------------
	Information on NAS Parallel Benchmarks is available at:
	http://www.nas.nasa.gov/Software/NPB/
	Authors: E. Barszcz
		P. Frederickson
		A. Woo
		M. Yarrow
	DASH version:
	Nicco Mietzsch <nicco.mietzsch@campus.lmu.de>
	CPP and OpenMP version:
			Dalvan Griebler <dalvangriebler@gmail.com>
			Júnior Löff <loffjh@gmail.com>
--------------------------------------------------------------------*/
#include <libdash.h>

#include <vector>
#include <iostream>
#include "npb-CPP.hpp"

#include "globals.hpp"

/* parameters */
#define T_BENCH	1
#define T_INIT	2
#define T_STL	3

/* global variables */
/* common /grid/ */
static int is1, is2, is3, ie1, ie2, ie3;

/* functions prototypes */
static void setup(int *n1, int *n2, int *n3, int lt);
static void mg3P(std::vector<dash::NArray<double, 3> > &u, dash::NArray<double, 3> &v, std::vector<dash::NArray<double, 3> > &r, double a[4], double c[4], int n1, int n2, int n3, int k);
static void psinv(dash::NArray<double, 3> &r, dash::NArray<double, 3> &u, int n1, int n2, int n3, double c[4], int k);
static void resid(dash::NArray<double, 3> &u, dash::NArray<double, 3> &v, dash::NArray<double, 3> &r, int n1, int n2, int n3, double a[4], int k);
static void rprj3(dash::NArray<double, 3> &r, int m1k, int m2k, int m3k, dash::NArray<double, 3> &s, int m1j, int m2j, int m3j, int k);
static void interp(dash::NArray<double, 3> &z, int mm1, int mm2, int mm3, dash::NArray<double, 3> &u, int n1, int n2, int n3, int k);
static void norm2u3(dash::NArray<double, 3> &r, int n1, int n2, int n3, double *rnm2, double *rnmu, int nx, int ny, int nz);
static void rep_nrm(dash::NArray<double, 3> &u, int n1, int n2, int n3, char *title, int kk);
static void comm3(dash::NArray<double, 3> &u, int n1, int n2, int n3);
static void zran3(dash::NArray<double, 3> &z, int n1, int n2, int n3, int nx, int ny, int k);
static void showall(dash::NArray<double, 3> &z, int n1, int n2, int n3);
static double power(double a, int n);
static void bubble(double ten[M][2], int j1[M][2], int j2[M][2], int j3[M][2], int m, int ind);
static void zero3(dash::NArray<double, 3> &z, int n1, int n2, int n3);
/*static void nonzero(double ***z, int n1, int n2, int n3);*/

/*--------------------------------------------------------------------
	program mg
c-------------------------------------------------------------------*/

int main(int argc, char *argv[]) {

	dash::init(&argc, &argv);
	/*-------------------------------------------------------------------------
	c k is the current level. It is passed down through subroutine args
	c and is NOT global. it is the current iteration
	c------------------------------------------------------------------------*/

	int k, it;
	double t, tinit, mflops;
	int nthreads = 1;

	/*-------------------------------------------------------------------------
	c These arrays are in common because they are quite large
	c and probably shouldn't be allocated on the stack. They
	c are always passed as subroutine args.
	c------------------------------------------------------------------------*/

	auto distspec = dash::DistributionSpec<3>( dash::BLOCKED, dash::NONE , dash::NONE);
	std::vector<dash::NArray<double, 3> > u;
	dash::NArray<double, 3> v;
	std::vector<dash::NArray<double, 3> > r;

	double a[4], c[4];

	double rnm2, rnmu;
	double epsilon = 1.0e-8;
	int n1, n2, n3, nit;
	double verify_value;
	boolean verified;

	int i, j, l;
	FILE *fp;

	timer_clear(T_BENCH);
	timer_clear(T_INIT);
	timer_clear(T_STL);

	nthreads = dash::size();

	timer_start(T_INIT);

	/*----------------------------------------------------------------------
	c Read in and broadcast input data
	c---------------------------------------------------------------------*/
	if(dash::myid() == 0) {
		printf("\n\n NAS Parallel Benchmarks 4.0 C++ DASH version" " - MG Benchmark\n\n");
		printf("\n\n Developed by: Dalvan Griebler <dalvan.griebler@acad.pucrs.br>\n");
		printf("\n\n DASH version by: Nicco Mietzsch <nicco.mietzsch@campus.lmu.de>\n");
	}

	fp = fopen("mg.input", "r");
	if (fp != NULL) {
		printf(" Reading from input file mg.input\n");
		if (fscanf(fp, "%d", &lt) != 1){
			printf(" Error in reading elements\n");
			exit(1);
		}
		while(fgetc(fp) != '\n');
		if (fscanf(fp, "%d%d%d", &nx[lt], &ny[lt], &nz[lt]) != 3){
			printf(" Error in reading elements\n");
			exit(1);
		}
		while(fgetc(fp) != '\n');
		if (fscanf(fp, "%d", &nit) != 1){
			printf(" Error in reading elements\n");
			exit(1);
		}
		while(fgetc(fp) != '\n');
		for (i = 0; i <= 7; i++) {
			if (fscanf(fp, "%d", &debug_vec[i]) != 1){
				printf(" Error in reading elements\n");
				exit(1);
			}
		}
		fclose(fp);
	} else {
		if(dash::myid() == 0) printf(" No input file. Using compiled defaults\n");

		lt = LT_DEFAULT;
		nit = NIT_DEFAULT;
		nx[lt] = NX_DEFAULT;
		ny[lt] = NY_DEFAULT;
		nz[lt] = NZ_DEFAULT;

		for (i = 0; i <= 7; i++) {
			debug_vec[i] = DEBUG_DEFAULT;
		}
	}

	if ( (nx[lt] != ny[lt]) || (nx[lt] != nz[lt]) ) {
		class_npb = 'U';
	} else if( nx[lt] == 32 && nit == 4 ) {
		class_npb = 'S';
	} else if( nx[lt] == 64 && nit == 40 ) {
		class_npb = 'W';
	} else if( nx[lt] == 256 && nit == 20 ) {
		class_npb = 'B';
	} else if( nx[lt] == 512 && nit == 20 ) {
		class_npb = 'C';
	} else if( nx[lt] == 256 && nit == 4 ) {
		class_npb = 'A';
	} else {
		class_npb = 'U';
	}

	/*--------------------------------------------------------------------
	c  Use these for debug info:
	c---------------------------------------------------------------------
	c	 debug_vec(0) = 1 !=> report all norms
	c	 debug_vec(1) = 1 !=> some setup information
	c	 debug_vec(1) = 2 !=> more setup information
	c	 debug_vec(2) = k => at level k or below, show result of resid
	c	 debug_vec(3) = k => at level k or below, show result of psinv
	c	 debug_vec(4) = k => at level k or below, show result of rprj
	c	 debug_vec(5) = k => at level k or below, show result of interp
	c	 debug_vec(6) = 1 => (unused)
	c	 debug_vec(7) = 1 => (unused)
	c-------------------------------------------------------------------*/

	a[0] = -8.0/3.0;
	a[1] =  0.0;
	a[2] =  1.0/6.0;
	a[3] =  1.0/12.0;

	if (class_npb == 'A' || class_npb == 'S' || class_npb =='W') {
		/*--------------------------------------------------------------------
		c	Coefficients for the S(a) smoother
		c-------------------------------------------------------------------*/
		c[0] =  -3.0/8.0;
		c[1] =  1.0/32.0;
		c[2] =  -1.0/64.0;
		c[3] =   0.0;
	} else {
		/*--------------------------------------------------------------------
		c	 Coefficients for the S(b) smoother
		c-------------------------------------------------------------------*/
		c[0] =  -3.0/17.0;
		c[1] =  1.0/33.0;
		c[2] =  -1.0/61.0;
		c[3] =   0.0;
	}

	lb = 1;

	setup(&n1,&n2,&n3,lt);

	u.resize(lt+1);
	for (l = lt; l >=1; l--) {
		u[l].allocate(m3[l],m2[l],m1[l], distspec);
	}

	v.allocate(m3[lt],m2[lt],m1[lt], distspec);

	r.resize(lt+1);
	for (l = lt; l >=1; l--) {
		r[l].allocate(m3[l],m2[l],m1[l], distspec);
	}

	zero3(u[lt],n1,n2,n3);
	zran3(v,n1,n2,n3,nx[lt],ny[lt],lt);

	norm2u3(v,n1,n2,n3,&rnm2,&rnmu,nx[lt],ny[lt],nz[lt]);

	if(dash::myid() == 0) {
		/*printf("\n norms of random v are\n");
		printf(" %4d%19.12e%19.12e\n", 0, rnm2, rnmu);
		printf(" about to evaluate resid, k= %d\n", lt);*/

		printf(" Size: %3dx%3dx%3d (class_npb %1c)\n", nx[lt], ny[lt], nz[lt], class_npb);
		printf(" Iterations: %3d\n", nit);
	}

	resid(u[lt],v,r[lt],n1,n2,n3,a,lt);
	norm2u3(r[lt],n1,n2,n3,&rnm2,&rnmu,nx[lt],ny[lt],nz[lt]);

	/*c---------------------------------------------------------------------
	c	 One iteration for startup
	c---------------------------------------------------------------------*/
	mg3P(u,v,r,a,c,n1,n2,n3,lt);
	resid(u[lt],v,r[lt],n1,n2,n3,a,lt);

	setup(&n1,&n2,&n3,lt);

	zero3(u[lt],n1,n2,n3);
	zran3(v,n1,n2,n3,nx[lt],ny[lt],lt);

	if(dash::myid() == 0) {
		timer_stop(T_INIT);

		timer_clear(T_STL);

		//std::clear();

		timer_start(T_BENCH);
	}

	resid(u[lt],v,r[lt],n1,n2,n3,a,lt);
	norm2u3(r[lt],n1,n2,n3,&rnm2,&rnmu,nx[lt],ny[lt],nz[lt]);

	for ( it = 1; it <= nit; it++) {
		mg3P(u,v,r,a,c,n1,n2,n3,lt);
		resid(u[lt],v,r[lt],n1,n2,n3,a,lt);
	}
	norm2u3(r[lt],n1,n2,n3,&rnm2,&rnmu,nx[lt],ny[lt],nz[lt]);

	if(dash::myid() == 0) {
		timer_stop(T_BENCH);

		t = timer_read(T_BENCH);
		tinit = timer_read(T_INIT);

		verified = FALSE;
		verify_value = 0.0;

		printf(" Initialization time: %15.3f seconds\n", tinit);
		printf(" Benchmark completed\n");

		if (class_npb != 'U') {
			if (class_npb == 'S') {
					verify_value = 0.530770700573e-04;
			} else if (class_npb == 'W') {
					verify_value = 0.250391406439e-17;  /* 40 iterations*/
				/*	0.183103168997d-044 iterations*/
			} else if (class_npb == 'A') {
					verify_value = 0.2433365309e-5;
				} else if (class_npb == 'B') {
					verify_value = 0.180056440132e-5;
				} else if (class_npb == 'C') {
					verify_value = 0.570674826298e-06;
			}

			if ( fabs( rnm2 - verify_value ) <= epsilon ) {
				verified = TRUE;
				printf(" VERIFICATION SUCCESSFUL\n");
				printf(" L2 Norm is %20.12e\n", rnm2);
				printf(" Error is   %20.12e\n", rnm2 - verify_value);
			} else {
				verified = FALSE;
				printf(" VERIFICATION FAILED\n");
				printf(" L2 Norm is			 %20.12e\n", rnm2);
				printf(" The correct L2 Norm is %20.12e\n", verify_value);
			}
		} else {
			verified = FALSE;
			printf(" Problem size unknown\n");
			printf(" NO VERIFICATION PERFORMED\n");
		}

		if ( t != 0.0 ) {
			int nn = nx[lt]*ny[lt]*nz[lt];
			mflops = 58.*nit*nn*1.0e-6 / t;
		} else {
			mflops = 0.0;
		}

		c_print_results((char*)"MG", class_npb, nx[lt], ny[lt], nz[lt], nit, nthreads, t, mflops, (char*)"		  floating point",
				verified, (char*)NPBVERSION, (char*)COMPILETIME, (char*)CS1, (char*)CS2, (char*)CS3, (char*)CS4, (char*)CS5, (char*)CS6, (char*)CS7);
		if(TIMERS_ENABLED == TRUE) printf(" time spent in STL: %15.3f seconds\n", timer_read(T_STL));

		//printf("\n mystl statistics:\n");
		//std::dump();
	}

	dash::finalize();

	return 0;
}

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void setup(int *n1, int *n2, int *n3, int lt) {

	/*--------------------------------------------------------------------
	c-------------------------------------------------------------------*/

	int k;

	for ( k = lt-1; k >= 1; k--) {
		nx[k] = nx[k+1]/2;
		ny[k] = ny[k+1]/2;
		nz[k] = nz[k+1]/2;
	}

	for (k = 1; k <= lt; k++) {
		m1[k] = nx[k]+2;
		m2[k] = nz[k]+2;
		m3[k] = ny[k]+2;
	}

	is1 = 1;
	ie1 = nx[lt];
	*n1 = nx[lt]+2;
	is2 = 1;
	ie2 = ny[lt];
	*n2 = ny[lt]+2;
	is3 = 1;
	ie3 = nz[lt];
	*n3 = nz[lt]+2;

	if (debug_vec[1] >=  1 && dash::myid() == 0) {
		printf(" in setup, \n");
		printf("  lt  nx  ny  nz  n1  n2  n3 is1 is2 is3 ie1 ie2 ie3\n");
		printf("%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d\n",lt,nx[lt],ny[lt],nz[lt],*n1,*n2,*n3,is1,is2,is3,ie1,ie2,ie3);
	}
}

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void mg3P(std::vector<dash::NArray<double, 3> > &u, dash::NArray<double, 3> &v, std::vector<dash::NArray<double, 3> > &r, double a[4], double c[4], int n1, int n2, int n3, int k) {

	/*--------------------------------------------------------------------
	c-------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
	c multigrid V-cycle routine
	c-------------------------------------------------------------------*/

	int j;

	/*--------------------------------------------------------------------
	c	 down cycle.
	c	 restrict the residual from the find grid to the coarse
	c-------------------------------------------------------------------*/

	for (k = lt; k >= lb+1; k--) {
		j = k-1;
		rprj3(r[k], m1[k], m2[k], m3[k], r[j], m1[j], m2[j], m3[j], k);
	}

	k = lb;
	/*--------------------------------------------------------------------
	c	 compute an approximate solution on the coarsest grid
	c-------------------------------------------------------------------*/
	zero3(u[k], m1[k], m2[k], m3[k]);
	psinv(r[k], u[k], m1[k], m2[k], m3[k], c, k);

	for (k = lb+1; k <= lt-1; k++) {
		j = k-1;
		/*--------------------------------------------------------------------
		c prolongate from level k-1  to k
		c-------------------------------------------------------------------*/
		zero3(u[k], m1[k], m2[k], m3[k]);
		interp(u[j], m1[j], m2[j], m3[j], u[k], m1[k], m2[k], m3[k], k);
		/*--------------------------------------------------------------------
		c compute residual for level k
		c-------------------------------------------------------------------*/
		resid(u[k], r[k], r[k], m1[k], m2[k], m3[k], a, k);
		/*--------------------------------------------------------------------
		c apply smoother
		c-------------------------------------------------------------------*/
		psinv(r[k], u[k], m1[k], m2[k], m3[k], c, k);
	}

	j = lt - 1;
	k = lt;
	interp(u[j], m1[j], m2[j], m3[j], u[lt], n1, n2, n3, k);
	resid(u[lt], v, r[lt], n1, n2, n3, a, k);
	psinv(r[lt], u[lt], n1, n2, n3, c, k);
}

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void psinv( dash::NArray<double, 3> &r, dash::NArray<double, 3> &u, int n1, int n2, int n3, double c[4], int k) {
	/*--------------------------------------------------------------------
	c-------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
	c	 psinv applies an approximate inverse as smoother:  u = u + Cr
	c
	c	 This  implementation costs  15A + 4M per result, where
	c	 A and M denote the costs of Addition and Multiplication.
	c	 Presuming coefficient c(3) is zero (the NPB assumes this,
	c	 but it is thus not a general case), 2A + 1M may be eliminated,
	c	 resulting in 13A + 3M.
	c	 Note that this vectorizes, and is also fine for cache
	c	 based machines.
	c-------------------------------------------------------------------*/
	double r1[n1], r2[n1];

	int myid = dash::myid();

	auto pattern = u.pattern();
	auto local_beg_gidx = pattern.coords(pattern.global(0));
  auto local_end_gidx = pattern.coords(pattern.global(pattern.local_size()-1));

	int topcoord = local_beg_gidx[0]-1;
	int bottomcoord = local_end_gidx[0]+1;
	int z_ext = u.local.extent(0);

	// double topplane[n2][n1];
	// double bottomplane[n2][n1];
	std::vector<double> topplane(n2*n1);
	std::vector<double> bottomplane(n2*n1);

	dash::Future<double*> fut_top;
	dash::Future<double*> fut_bot;

	if(topcoord >= 0 && z_ext > 0) {
		fut_top = dash::copy_async(r.begin()+n2*n1*topcoord, r.begin()+n2*n1*(topcoord+1), &topplane[0]);
	}

	if(bottomcoord < n3 && z_ext > 0) {
		fut_bot = dash::copy_async(r.begin()+n2*n1*bottomcoord, r.begin()+n2*n1*(bottomcoord+1), &bottomplane[0]);
	}

	for(int i3 = 1; i3 < z_ext-1; i3++) {
		for (int i2 = 1; i2 < n2-1; i2++) {
			for (int i1 = 0; i1 < n1; i1++) {
				r1[i1] = r.local(i3,i2-1,i1) + r.local(i3,i2+1,i1) + r.local(i3-1,i2,i1) + r.local(i3+1,i2,i1);
				r2[i1] = r.local(i3-1,i2-1,i1) + r.local(i3-1,i2+1,i1) + r.local(i3+1,i2-1,i1) + r.local(i3+1,i2+1,i1);
			}
			for (int i1 = 1; i1 < n1-1; i1++) {
				u.local(i3,i2,i1) = u.local(i3,i2,i1)
				+ c[0] * r.local(i3,i2,i1)
				+ c[1] * ( r.local(i3,i2,i1-1) + r.local(i3,i2,i1+1) + r1[i1] )
				+ c[2] * ( r2[i1] + r1[i1-1] + r1[i1+1] );
				//--------------------------------------------------------------------
				//c  Assume c(3) = 0	(Enable line below if c(3) not= 0)
				//c-------------------------------------------------------------------
				//c > + c(3) * ( r2(i1-1) + r2(i1+1) )
				//c-------------------------------------------------------------------
			}
		}
	}

	if(z_ext > 1) {
		if(topcoord >= 0) {
			int i3 = 0;
			fut_top.wait();

			for (int i2 = 1; i2 < n2-1; i2++) {
				for (int i1 = 0; i1 < n1; i1++) {
					r1[i1] = r.local(i3,i2-1,i1) + r.local(i3,i2+1,i1) + topplane[i2*n2+i1] + r.local(i3+1,i2,i1);
					r2[i1] = topplane[(i2-1)*n2+i1] + topplane[(i2+1)*n2+i1] + r.local(i3+1,i2-1,i1) + r.local(i3+1,i2+1,i1);
				}
				for (int i1 = 1; i1 < n1-1; i1++) {
					u.local(i3,i2,i1) = u.local(i3,i2,i1)
					+ c[0] * r.local(i3,i2,i1)
					+ c[1] * ( r.local(i3,i2,i1-1) + r.local(i3,i2,i1+1) + r1[i1] )
					+ c[2] * ( r2[i1] + r1[i1-1] + r1[i1+1] );
					//--------------------------------------------------------------------
					//c  Assume c(3) = 0	(Enable line below if c(3) not= 0)
					//c-------------------------------------------------------------------
					//c > + c(3) * ( r2(i1-1) + r2(i1+1) )
					//c-------------------------------------------------------------------
				}
			}
		}
		if(bottomcoord < n3) {
			int i3 = z_ext-1;
			fut_bot.wait();

			for (int i2 = 1; i2 < n2-1; i2++) {
				for (int i1 = 0; i1 < n1; i1++) {
					r1[i1] = r.local(i3,i2-1,i1) + r.local(i3,i2+1,i1) + r.local(i3-1,i2,i1) + bottomplane[i2*n2+i1];
					r2[i1] = r.local(i3-1,i2-1,i1) + r.local(i3-1,i2+1,i1) + bottomplane[(i2-1)*n2+i1] + bottomplane[(i2+1)*n2+i1];
				}
				for (int i1 = 1; i1 < n1-1; i1++) {
					u.local(i3,i2,i1) = u.local(i3,i2,i1)
					+ c[0] * r.local(i3,i2,i1)
					+ c[1] * ( r.local(i3,i2,i1-1) + r.local(i3,i2,i1+1) + r1[i1] )
					+ c[2] * ( r2[i1] + r1[i1-1] + r1[i1+1] );
					//--------------------------------------------------------------------
					//c  Assume c(3) = 0	(Enable line below if c(3) not= 0)
					//c-------------------------------------------------------------------
					//c > + c(3) * ( r2(i1-1) + r2(i1+1) )
					//c-------------------------------------------------------------------
				}
			}
		}
	} else {
		if(0 <= topcoord && bottomcoord < n3) {
			int i3 = 0;
			fut_top.wait();
			fut_bot.wait();

			for (int i2 = 1; i2 < n2-1; i2++) {
				for (int i1 = 0; i1 < n1; i1++) {
					r1[i1] = r.local(i3,i2-1,i1) + r.local(i3,i2+1,i1) + topplane[i2*n2+i1] + bottomplane[i2*n2+i1];
					r2[i1] = topplane[(i2-1)*n2+i1] + topplane[(i2+1)*n2+i1] + bottomplane[(i2-1)*n2+i1] + bottomplane[(i2+1)*n2+i1];
				}
				for (int i1 = 1; i1 < n1-1; i1++) {
					u.local(i3,i2,i1) = u.local(i3,i2,i1)
					+ c[0] * r.local(i3,i2,i1)
					+ c[1] * ( r.local(i3,i2,i1-1) + r.local(i3,i2,i1+1) + r1[i1] )
					+ c[2] * ( r2[i1] + r1[i1-1] + r1[i1+1] );
					//--------------------------------------------------------------------
					//c  Assume c(3) = 0	(Enable line below if c(3) not= 0)
					//c-------------------------------------------------------------------
					//c > + c(3) * ( r2(i1-1) + r2(i1+1) )
					//c-------------------------------------------------------------------
				}
			}
		}
	}
		/*--------------------------------------------------------------------
		c	 exchange boundary points
		c-------------------------------------------------------------------*/
		comm3(u,n1,n2,n3);

		if (debug_vec[0] >= 1 ) {
			rep_nrm(u,n1,n2,n3,(char*)"   psinv",k);
		}

		if ( debug_vec[3] >= k ) {
			showall(u,n1,n2,n3);
		}
	}

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void resid( dash::NArray<double, 3> &u, dash::NArray<double, 3> &v, dash::NArray<double, 3> &r, int n1, int n2, int n3, double a[4], int k ) {
	/*--------------------------------------------------------------------
	c-------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
	c	 resid computes the residual:  r = v - Au
	c
	c	 This  implementation costs  15A + 4M per result, where
	c	 A and M denote the costs of Addition (or Subtraction) and
	c	 Multiplication, respectively.
	c	 Presuming coefficient a(1) is zero (the NPB assumes this,
	c	 but it is thus not a general case), 3A + 1M may be eliminated,
	c	 resulting in 12A + 3M.
	c	 Note that this vectorizes, and is also fine for cache
	c	 based machines.
	c-------------------------------------------------------------------*/
	// double u1[n1], u2[n1];
	std::vector<double> u1(n1);
	std::vector<double> u2(n1);

	int myid = dash::myid();

	auto pattern = r.pattern();
	auto local_beg_gidx = pattern.coords(pattern.global(0));
	auto local_end_gidx = pattern.coords(pattern.global(pattern.local_size()-1));

	int topcoord = local_beg_gidx[0]-1;
	int bottomcoord = local_end_gidx[0]+1;
	int z_ext = u.local.extent(0);

	// double topplane[n2][n1];
	// double bottomplane[n2][n1];
	std::vector<double> topplane(n2*n1);
	std::vector<double> bottomplane(n2*n1);
	// printf("init done.\n");

	dash::Future<double*> fut_top;
	dash::Future<double*> fut_bot;

	if(topcoord >= 0 && z_ext > 0) {
		fut_top = dash::copy_async(u.begin()+n2*n1*topcoord, u.begin()+n2*n1*(topcoord+1), &topplane[0]);
	}

	if(bottomcoord < n3 && z_ext > 0) {
		fut_bot = dash::copy_async(u.begin()+n2*n1*bottomcoord, u.begin()+n2*n1*(bottomcoord+1), &bottomplane[0]);
	}
	// printf("copy done.\n");

	for(int i3 = 1; i3 < z_ext-1; i3++) {
		// printf("Unit %d trying local plane %d of %d\n", (int) dash::myid(), i3, z_ext);
		for (int i2 = 1; i2 < n2-1; i2++) {
			for (int i1 = 0; i1 < n1; i1++) {
				u1[i1] = u.local(i3,i2-1,i1) + u.local(i3,i2+1,i1) + u.local(i3-1,i2,i1) + u.local(i3+1,i2,i1);
				u2[i1] = u.local(i3-1,i2-1,i1) + u.local(i3-1,i2+1,i1) + u.local(i3+1,i2-1,i1) + u.local(i3+1,i2+1,i1);
			}
			for (int i1 = 1; i1 < n1-1; i1++) {
				r.local(i3,i2,i1) = v.local(i3,i2,i1)
				 - a[0] * u.local(i3,i2,i1)
				//--------------------------------------------------------------------
				//c  Assume a(1) = 0	  (Enable 2 lines below if a(1) not= 0)
				//c-------------------------------------------------------------------
				//c > - a(1) * ( u(i1-1,i2,i3) + u(i1+1,i2,i3)
				//c > + u1(i1) )
				//c-------------------------------------------------------------------
				 - a[2] * ( u2[i1] + u1[i1-1] + u1[i1+1] )
				 - a[3] * ( u2[i1-1] + u2[i1+1] );
			}
		}
	}
	// printf("Main iter done.\n");

	if(z_ext > 1) {
		if(topcoord >= 0) {
			int i3 = 0;
			fut_top.wait();

			for (int i2 = 1; i2 < n2-1; i2++) {
				for (int i1 = 0; i1 < n1; i1++) {
					u1[i1] = u.local(i3,i2-1,i1) + u.local(i3,i2+1,i1) + topplane[i2*n2+i1] + u.local(i3+1,i2,i1);
					u2[i1] = topplane[(i2-1)*n2+i1] + topplane[(i2+1)*n2+i1] + u.local(i3+1,i2-1,i1) + u.local(i3+1,i2+1,i1);
				}
				for (int i1 = 1; i1 < n1-1; i1++) {
					r.local(i3,i2,i1) = v.local(i3,i2,i1)
					 - a[0] * u.local(i3,i2,i1)
					//--------------------------------------------------------------------
					//c  Assume a(1) = 0	  (Enable 2 lines below if a(1) not= 0)
					//c-------------------------------------------------------------------
					//c > - a(1) * ( u(i1-1,i2,i3) + u(i1+1,i2,i3)
					//c > + u1(i1) )
					//c-------------------------------------------------------------------
					 - a[2] * ( u2[i1] + u1[i1-1] + u1[i1+1] )
					 - a[3] * ( u2[i1-1] + u2[i1+1] );
				}
			}
		}
		if(bottomcoord < n3) {
			int i3 = z_ext-1;
			fut_bot.wait();

			for (int i2 = 1; i2 < n2-1; i2++) {
				for (int i1 = 0; i1 < n1; i1++) {
					u1[i1] = u.local(i3,i2-1,i1) + u.local(i3,i2+1,i1) + u.local(i3-1,i2,i1) + bottomplane[i2*n2+i1];
					u2[i1] = u.local(i3-1,i2-1,i1) + u.local(i3-1,i2+1,i1) + bottomplane[(i2-1)*n2+i1] + bottomplane[(i2+1)*n2+i1];
				}
				for (int i1 = 1; i1 < n1-1; i1++) {
					r.local(i3,i2,i1) = v.local(i3,i2,i1)
					 - a[0] * u.local(i3,i2,i1)
					//--------------------------------------------------------------------
					//c  Assume a(1) = 0	  (Enable 2 lines below if a(1) not= 0)
					//c-------------------------------------------------------------------
					//c > - a(1) * ( u(i1-1,i2,i3) + u(i1+1,i2,i3)
					//c > + u1(i1) )
					//c-------------------------------------------------------------------
					 - a[2] * ( u2[i1] + u1[i1-1] + u1[i1+1] )
					 - a[3] * ( u2[i1-1] + u2[i1+1] );
				}
			}
		}
	} else {
		if(0 <= topcoord && bottomcoord < n3) {
			int i3 = 0;
			fut_top.wait();
			fut_bot.wait();

			for (int i2 = 1; i2 < n2-1; i2++) {
				for (int i1 = 0; i1 < n1; i1++) {
					u1[i1] = u.local(i3,i2-1,i1) + u.local(i3,i2+1,i1) + topplane[i2*n2+i1] + bottomplane[i2*n2+i1];
					u2[i1] = topplane[(i2-1)*n2+i1] + topplane[(i2+1)*n2+i1] + bottomplane[(i2-1)*n2+i1] + bottomplane[(i2+1)*n2+i1];
				}
				for (int i1 = 1; i1 < n1-1; i1++) {
					r.local(i3,i2,i1) = v.local(i3,i2,i1)
					 - a[0] * u.local(i3,i2,i1)
					//--------------------------------------------------------------------
					//c  Assume a(1) = 0	  (Enable 2 lines below if a(1) not= 0)
					//c-------------------------------------------------------------------
					//c > - a(1) * ( u(i1-1,i2,i3) + u(i1+1,i2,i3)
					//c > + u1(i1) )
					//c-------------------------------------------------------------------
					 - a[2] * ( u2[i1] + u1[i1-1] + u1[i1+1] )
					 - a[3] * ( u2[i1-1] + u2[i1+1] );
				}
			}
		}
	}

	/*--------------------------------------------------------------------
	c	 exchange boundary data
	c--------------------------------------------------------------------*/
	comm3(r,n1,n2,n3);

	if (debug_vec[0] >= 1 ) {
		rep_nrm(r,n1,n2,n3,(char*)"   resid",k);
	}

	if ( debug_vec[2] >= k ) {
		showall(r,n1,n2,n3);
	}
}

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void rprj3( dash::NArray<double, 3> &r, int m1k, int m2k, int m3k, dash::NArray<double, 3> &s, int m1j, int m2j, int m3j, int k ) {
	/*--------------------------------------------------------------------
		c-------------------------------------------------------------------*/

		/*--------------------------------------------------------------------
		c	 rprj3 projects onto the next coarser grid,
		c	 using a trilinear Finite Element projection:  s = r' = P r
		c
		c	 This  implementation costs  20A + 4M per result, where
		c	 A and M denote the costs of Addition and Multiplication.
		c	 Note that this vectorizes, and is also fine for cache
		c	 based machines.
		c-------------------------------------------------------------------*/

		int d1, d2, d3;

		if (m1k == 3) {
			d1 = 2;
		} else {
			d1 = 1;
		}

		if (m2k == 3) {
			d2 = 2;
		} else {
			d2 = 1;
		}

		if (m3k == 3) {
			d3 = 2;
		} else {
			d3 = 1;
		}

		int j2, j1, i3, i2, i1;
		double x1[m1k], y1[m1k], x2, y2;

		auto s_pattern = s.pattern();
		auto s_local_beg_gidx = s_pattern.coords(s_pattern.global(0));
	  auto s_local_end_gidx = s_pattern.coords(s_pattern.global(s_pattern.local_size()-1));

		auto r_pattern = r.pattern();
		auto r_local_beg_gidx = r_pattern.coords(r_pattern.global(0));
	  auto r_local_end_gidx = r_pattern.coords(r_pattern.global(r_pattern.local_size()-1));
		int r_offset = r_local_beg_gidx[0];
		int r_x = (int) r.extent(2);
		int r_y = (int) r.extent(1);

		int start = 0;
		if(s_local_beg_gidx[0] == 0) start++;

		int end = s.local.extent(0);
		if(s_local_end_gidx[0] == s.extent(0)-1) end--;

		int r_planes = 2*(end-start)+1;

		std::vector<std::vector<double> > r_local_v(r_planes);
    for(int i = 0; i < r_planes; i++) {
      r_local_v[i] = std::vector<double> (r_y*r_x);
    }
		std::vector<double*> r_local(r_planes);
		std::vector<bool> is_local(r_planes);
		int r_idx = 0;
		int r_psize = r_y*r_x;

		dash::Future<double*> futs[r_planes];

		for(int j3l = start; j3l < end; j3l++){
			int j3 = s_local_beg_gidx[0]+j3l;
			i3 = 2*j3-d3;
			if(r(i3,0,0).is_local()) {
				r_local[r_idx] = r.lbegin()+r_psize*(i3-r_offset);
				is_local[r_idx] = true;
				futs[r_idx] = dash::copy_async(r.begin()+r_psize*i3, r.begin()+r_psize*i3, r_local_v[r_idx].data());
			} else {
				futs[r_idx] = dash::copy_async(r.begin()+r_psize*i3, r.begin()+r_psize*(i3+1), r_local_v[r_idx].data());
				r_local[r_idx] = r_local_v[r_idx].data();
			}
			r_idx++;
			if(r(i3+1,0,0).is_local()) {
				r_local[r_idx] = r.lbegin()+r_psize*(i3+1-r_offset);
				is_local[r_idx] = true;
				futs[r_idx] = dash::copy_async(r.begin()+r_psize*(i3+1), r.begin()+r_psize*(i3+1), r_local_v[r_idx].data());
			} else {
				futs[r_idx] = dash::copy_async(r.begin()+r_psize*(i3+1), r.begin()+r_psize*(i3+2), r_local_v[r_idx].data());
				r_local[r_idx] = r_local_v[r_idx].data();
			}
			r_idx++;
		}


		if(start < end) {
			i3 = 2*(s_local_beg_gidx[0]+end)-d3;
			if(r(i3,0,0).is_local()) {
				r_local[r_idx] = r.lbegin()+r_psize*(i3-r_offset);
				is_local[r_idx] = true;
				futs[r_idx] = dash::copy_async(r.begin()+r_psize*i3, r.begin()+r_psize*i3, r_local_v[r_idx].data());
			} else {
				futs[r_idx] = dash::copy_async(r.begin()+r_psize*i3, r.begin()+r_psize*(i3+1), r_local_v[r_idx].data());
				r_local[r_idx] = r_local_v[r_idx].data();
			}
		}

		//First, do the local calculations...
		r_idx = 0;

		for(int j3l = start; j3l < end; j3l++) {
			if(is_local[r_idx+0] && is_local[r_idx+1] && is_local[r_idx+2]) {
				int j3 = s_local_beg_gidx[0]+j3l;
				i3 = 2*j3-d3;
				//C	i3 = 2*j3-1

				for (j2 = 1; j2 < m2j-1; j2++) {
					i2 = 2*j2-d2;
					//C  i2 = 2*j2-1

					for (j1 = 1; j1 < m1j; j1++) {
						i1 = 2*j1-d1;
					//C	i1 = 2*j1-1
					  x1[i1] = r_local[r_idx+1][i2*r_x+i1] + r_local[r_idx+1][(i2+2)*r_x+i1] + r_local[r_idx+0][(i2+1)*r_x+i1] + r_local[r_idx+2][(i2+1)*r_x+i1];
						y1[i1] = r_local[r_idx+0][i2*r_x+i1] + r_local[r_idx+2][i2*r_x+i1] + r_local[r_idx+0][(i2+2)*r_x+i1] + r_local[r_idx+2][(i2+2)*r_x+i1];
					}

					for (j1 = 1; j1 < m1j-1; j1++) {
						i1 = 2*j1-d1;
					//C	i1 = 2*j1-1
						y2 = r_local[r_idx+0][i2*r_x+i1+1] + r_local[r_idx+2][i2*r_x+i1+1] + r_local[r_idx+0][(i2+2)*r_x+i1+1] + r_local[r_idx+2][(i2+2)*r_x+i1+1];
						x2 = r_local[r_idx+1][i2*r_x+i1+1] + r_local[r_idx+1][(i2+2)*r_x+i1+1] + r_local[r_idx+0][(i2+1)*r_x+i1+1] + r_local[r_idx+2][(i2+1)*r_x+i1+1];
						s.local(j3l,j2,j1) =
							0.5 * r_local[r_idx+1][(i2+1)*r_x+i1+1]
							+ 0.25 * ( r_local[r_idx+1][(i2+1)*r_x+i1] + r_local[r_idx+1][(i2+1)*r_x+i1+2] + x2)
							+ 0.125 * ( x1[i1] + x1[i1+2] + y2)
							+ 0.0625 * ( y1[i1] + y1[i1+2] );
					}
				}
			}
			r_idx = r_idx + 2;
		}

		//...then do the non-local calculations.
		r_idx = 0;

		for(int j3l = start; j3l < end; j3l++) {
			if(!(is_local[r_idx+0] && is_local[r_idx+1] && is_local[r_idx+2])) {
				int j3 = s_local_beg_gidx[0]+j3l;
				i3 = 2*j3-d3;
				//C	i3 = 2*j3-1

				futs[r_idx+0].wait();
				futs[r_idx+1].wait();
				futs[r_idx+2].wait();

				for (j2 = 1; j2 < m2j-1; j2++) {
					i2 = 2*j2-d2;
					//C  i2 = 2*j2-1

					for (j1 = 1; j1 < m1j; j1++) {
						i1 = 2*j1-d1;
					//C	i1 = 2*j1-1
					  x1[i1] = r_local[r_idx+1][i2*r_x+i1] + r_local[r_idx+1][(i2+2)*r_x+i1] + r_local[r_idx+0][(i2+1)*r_x+i1] + r_local[r_idx+2][(i2+1)*r_x+i1];
						y1[i1] = r_local[r_idx+0][i2*r_x+i1] + r_local[r_idx+2][i2*r_x+i1] + r_local[r_idx+0][(i2+2)*r_x+i1] + r_local[r_idx+2][(i2+2)*r_x+i1];
					}

					for (j1 = 1; j1 < m1j-1; j1++) {
						i1 = 2*j1-d1;
					//C	i1 = 2*j1-1
						y2 = r_local[r_idx+0][i2*r_x+i1+1] + r_local[r_idx+2][i2*r_x+i1+1] + r_local[r_idx+0][(i2+2)*r_x+i1+1] + r_local[r_idx+2][(i2+2)*r_x+i1+1];
						x2 = r_local[r_idx+1][i2*r_x+i1+1] + r_local[r_idx+1][(i2+2)*r_x+i1+1] + r_local[r_idx+0][(i2+1)*r_x+i1+1] + r_local[r_idx+2][(i2+1)*r_x+i1+1];
						s.local(j3l,j2,j1) =
							0.5 * r_local[r_idx+1][(i2+1)*r_x+i1+1]
							+ 0.25 * ( r_local[r_idx+1][(i2+1)*r_x+i1] + r_local[r_idx+1][(i2+1)*r_x+i1+2] + x2)
							+ 0.125 * ( x1[i1] + x1[i1+2] + y2)
							+ 0.0625 * ( y1[i1] + y1[i1+2] );
					}
				}
			}
			r_idx = r_idx + 2;
		}

		comm3(s,m1j,m2j,m3j);

		if (debug_vec[0] >= 1 ) {
			rep_nrm(s,m1j,m2j,m3j,(char*)"   rprj3",k-1);
		}

		if (debug_vec[4] >= k ) {
			showall(s,m1j,m2j,m3j);
		}
	}

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void interp( dash::NArray<double, 3> &z, int mm1, int mm2, int mm3, dash::NArray<double, 3> &u, int n1, int n2, int n3, int k ) {
	/*--------------------------------------------------------------------
	c-------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
	c	 interp adds the trilinear interpolation of the correction
	c	 from the coarser grid to the current approximation:  u = u + Qu'
	c
	c	 Observe that this  implementation costs  16A + 4M, where
	c	 A and M denote the costs of Addition and Multiplication.
	c	 Note that this vectorizes, and is also fine for cache
	c	 based machines.  Vector machines may get slightly better
	c	 performance however, with 8 separate "do i1" loops, rather than 4.
	c-------------------------------------------------------------------*/

	int d1, d2, d3, t1, t2, t3;

	/*
	c note that m = 1037 in globals.h but for this only need to be
	c 535 to handle up to 1024^3
	c integer m
	c parameter( m=535 )
	*/

	if ( n1 != 3 && n2 != 3 && n3 != 3 ) {

		double z1[mm1], z2[mm1], z3[mm1];

		auto pattern = u.pattern();
		auto local_beg_gidx = pattern.coords(pattern.global(0));
	  auto local_end_gidx = pattern.coords(pattern.global(pattern.local_size()-1));

		auto z_pattern = z.pattern();
		auto z_local_beg_gidx = z_pattern.coords(z_pattern.global(0));
	  auto z_local_end_gidx = z_pattern.coords(z_pattern.global(z_pattern.local_size()-1));
		int z_offset = z_local_beg_gidx[0];

		int u_offset = local_beg_gidx[0];

		int start = 0 - u_offset%2;
		int u_planes = u.local.extent(0);
		int z_size = mm2*mm1;

		if((int) u_planes > 0) {
			int z_local_size = u_planes-start+((u_planes+start)%2);

			std::vector<std::vector<double> > z_local_v(z_local_size);
	    for(int i = 0; i < z_local_size; i++) {
	      z_local_v[i] = std::vector<double> (mm2*mm1);
	    }
			std::vector<double*> z_local(z_local_size);
			std::vector<bool> is_local(z_local_size);

			dash::Future<double*> futs[z_local_size];

			for(int j3 = start; j3 < u_planes; j3+=2) {
				int i3 = (u_offset + j3)/2;
				if(z(i3,0,0).is_local()) {
					z_local[j3-start] = z.lbegin()+z_size*(i3-z_offset);
					is_local[j3-start] = true;
					futs[j3-start] = dash::copy_async(z.begin()+z_size*i3, z.begin()+z_size*i3, z_local_v[j3-start].data());
				} else {
					futs[j3-start] = dash::copy_async(z.begin()+z_size*i3, z.begin()+z_size*(i3+1), z_local_v[j3-start].data());
					z_local[j3-start] = z_local_v[j3-start].data();
				}
				if(z(i3+1,0,0).is_local()) {
					z_local[j3-start+1] = z.lbegin()+z_size*(i3+1-z_offset);
					futs[j3-start+1] = dash::copy_async(z.begin()+z_size*(i3+1), z.begin()+z_size*(i3+1), z_local_v[j3-start+1].data());
					is_local[j3-start+1] = true;
				} else {
						futs[j3-start+1] = dash::copy_async(z.begin()+z_size*(i3+1), z.begin()+z_size*(i3+2), z_local_v[j3-start+1].data());
						z_local[j3-start+1] = z_local_v[j3-start+1].data();
				}
			}

			//First, do the local calculations...
			for(int j3 = start; j3 < u_planes; j3+=2) {
				// int i3 = (u_offset + j3)/2;
				if(is_local[j3-start] && is_local[j3-start+1]) {

					for (int i2 = 0; i2 < mm2-1; i2++) {
						for (int i1 = 0; i1 < mm1; i1++) {
							z1[i1] = z_local[j3-start][(i2+1)*mm2+i1] + z_local[j3-start][i2*mm2+i1];
							z2[i1] = z_local[j3-start+1][i2*mm2+i1] + z_local[j3-start][i2*mm2+i1];
							z3[i1] = z_local[j3-start+1][(i2+1)*mm2+i1] + z_local[j3-start+1][i2*mm2+i1] + z1[i1];
						}
						if(j3 >= 0) {
							for (int i1 = 0; i1 < mm1-1; i1++) {
								u.local(j3,2*i2,2*i1)     = u.local(j3,2*i2,2*i1) + z_local[j3-start][i2*mm2+i1];
								u.local(j3,2*i2,2*i1+1)   = u.local(j3,2*i2,2*i1+1) + 0.5*(z_local[j3-start][i2*mm2+i1+1]+z_local[j3-start][i2*mm2+i1]);
								u.local(j3,2*i2+1,2*i1)   = u.local(j3,2*i2+1,2*i1) + 0.5 * z1[i1];
								u.local(j3,2*i2+1,2*i1+1) = u.local(j3,2*i2+1,2*i1+1) + 0.25*( z1[i1] + z1[i1+1] );
							}
						}
						if(j3+1 < u_planes) {
							for (int i1 = 0; i1 < mm1-1; i1++) {
								u.local(j3+1,2*i2,2*i1)     = u.local(j3+1,2*i2,2*i1) + 0.5 * z2[i1];
								u.local(j3+1,2*i2,2*i1+1)   = u.local(j3+1,2*i2,2*i1+1) + 0.25*( z2[i1] + z2[i1+1] );
								u.local(j3+1,2*i2+1,2*i1)   = u.local(j3+1,2*i2+1,2*i1) + 0.25* z3[i1];
								u.local(j3+1,2*i2+1,2*i1+1) = u.local(j3+1,2*i2+1,2*i1+1) + 0.125*( z3[i1] + z3[i1+1] );
							}
						}
					}
				}
			}

			//...then do the non-local calculations.
			for(int j3 = start; j3 < u_planes; j3+=2) {
				if(!(is_local[j3-start] && is_local[j3-start+1])) {
					// int i3 = (u_offset + j3)/2;
					futs[j3-start].wait();
					futs[j3-start+1].wait();

					for (int i2 = 0; i2 < mm2-1; i2++) {
						for (int i1 = 0; i1 < mm1; i1++) {
							z1[i1] = z_local[j3-start][(i2+1)*mm2+i1] + z_local[j3-start][i2*mm2+i1];
							z2[i1] = z_local[j3-start+1][i2*mm2+i1] + z_local[j3-start][i2*mm2+i1];
							z3[i1] = z_local[j3-start+1][(i2+1)*mm2+i1] + z_local[j3-start+1][i2*mm2+i1] + z1[i1];
						}
						if(j3 >= 0) {
							for (int i1 = 0; i1 < mm1-1; i1++) {
								u.local(j3,2*i2,2*i1)     = u.local(j3,2*i2,2*i1) + z_local[j3-start][i2*mm2+i1];
								u.local(j3,2*i2,2*i1+1)   = u.local(j3,2*i2,2*i1+1) + 0.5*(z_local[j3-start][i2*mm2+i1+1]+z_local[j3-start][i2*mm2+i1]);
								u.local(j3,2*i2+1,2*i1)   = u.local(j3,2*i2+1,2*i1) + 0.5 * z1[i1];
								u.local(j3,2*i2+1,2*i1+1) = u.local(j3,2*i2+1,2*i1+1) + 0.25*( z1[i1] + z1[i1+1] );
							}
						}
						if(j3+1 < u_planes) {
							for (int i1 = 0; i1 < mm1-1; i1++) {
								u.local(j3+1,2*i2,2*i1)     = u.local(j3+1,2*i2,2*i1) + 0.5 * z2[i1];
								u.local(j3+1,2*i2,2*i1+1)   = u.local(j3+1,2*i2,2*i1+1) + 0.25*( z2[i1] + z2[i1+1] );
								u.local(j3+1,2*i2+1,2*i1)   = u.local(j3+1,2*i2+1,2*i1) + 0.25* z3[i1];
								u.local(j3+1,2*i2+1,2*i1+1) = u.local(j3+1,2*i2+1,2*i1+1) + 0.125*( z3[i1] + z3[i1+1] );
							}
						}
					}
				}
			}
		}

	} else { //this gets never called for our benchmark sizes, so we didn't parallelize it.

		if (n1 == 3) {
			d1 = 2;
			t1 = 1;
		} else {
			d1 = 1;
			t1 = 0;
		}
		if (n2 == 3) {
			d2 = 2;
			t2 = 1;
		} else {
			d2 = 1;
			t2 = 0;
		}
		if (n3 == 3) {
			d3 = 2;
			t3 = 1;
		} else {
			d3 = 1;
			t3 = 0;
		}
		if( 0 == dash::myid()) {

		for(int i3 = d3; i3 <= mm3-1; i3++) {
			for (int i2 = d2; i2 <= mm2-1; i2++) {
				for (int i1 = d1; i1 <= mm1-1; i1++) {
					u[2*i3-d3-1][2*i2-d2-1][2*i1-d1-1] =
					u[2*i3-d3-1][2*i2-d2-1][2*i1-d1-1]
					+z[i3-1][i2-1][i1-1];
				}
				for (int i1 = 1; i1 <= mm1-1; i1++) {
					u[2*i3-d3-1][2*i2-d2-1][2*i1-t1-1] =
					u[2*i3-d3-1][2*i2-d2-1][2*i1-t1-1]
					+0.5*(z[i3-1][i2-1][i1]+z[i3-1][i2-1][i1-1]);
				}
			}
			for (int i2 = 1; i2 <= mm2-1; i2++) {
				for (int i1 = d1; i1 <= mm1-1; i1++) {
					u[2*i3-d3-1][2*i2-t2-1][2*i1-d1-1] =
					u[2*i3-d3-1][2*i2-t2-1][2*i1-d1-1]
					+0.5*(z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
				}
				for (int i1 = 1; i1 <= mm1-1; i1++) {
					u[2*i3-d3-1][2*i2-t2-1][2*i1-t1-1] =
					u[2*i3-d3-1][2*i2-t2-1][2*i1-t1-1]
					+0.25*(z[i3-1][i2][i1]+z[i3-1][i2-1][i1] + z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
				}
			}
		}

		for(int i3 = 1; i3 <= mm3-1; i3++) {
			for (int i2 = d2; i2 <= mm2-1; i2++) {
				for (int i1 = d1; i1 <= mm1-1; i1++) {
					u[2*i3-t3-1][2*i2-d2-1][2*i1-d1-1] =
					u[2*i3-t3-1][2*i2-d2-1][2*i1-d1-1]
					+0.5*(z[i3][i2-1][i1-1]+z[i3-1][i2-1][i1-1]);
				}
				for (int i1 = 1; i1 <= mm1-1; i1++) {
					u[2*i3-t3-1][2*i2-d2-1][2*i1-t1-1] =
					u[2*i3-t3-1][2*i2-d2-1][2*i1-t1-1]
					+0.25*(z[i3][i2-1][i1]+z[i3][i2-1][i1-1] + z[i3-1][i2-1][i1]+z[i3-1][i2-1][i1-1]);
				}
			}
			for (int i2 = 1; i2 <= mm2-1; i2++) {
				for (int i1 = d1; i1 <= mm1-1; i1++) {
					u[2*i3-t3-1][2*i2-t2-1][2*i1-d1-1] =
					u[2*i3-t3-1][2*i2-t2-1][2*i1-d1-1]
					+0.25*(z[i3][i2][i1-1]+z[i3][i2-1][i1-1] + z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
				}
				for (int i1 = 1; i1 <= mm1-1; i1++) {
					u[2*i3-t3-1][2*i2-t2-1][2*i1-t1-1] =
					u[2*i3-t3-1][2*i2-t2-1][2*i1-t1-1]
					+0.125*(z[i3][i2][i1]+z[i3][i2-1][i1]
						+z[i3][i2][i1-1]+z[i3][i2-1][i1-1]
						+z[i3-1][i2][i1]+z[i3-1][i2-1][i1]
						+z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
				}
			}
		}


	{
		if (debug_vec[0] >= 1 ) {
			rep_nrm(z,mm1,mm2,mm3,(char*)"z: inter",k-1);
			rep_nrm(u,n1,n2,n3,(char*)"u: inter",k);
		}
		if ( debug_vec[5] >= k ) {
			showall(z,mm1,mm2,mm3);
			showall(u,n1,n2,n3);
		}
	}}
}
dash::barrier();
}

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void norm2u3( dash::NArray<double, 3> &r, int n1, int n2, int n3, double *rnm2, double *rnmu, int nx, int ny, int nz) {
	/*--------------------------------------------------------------------
	c-------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
	c	 norm2u3 evaluates approximations to the L2 norm and the
	c	 uniform (or L-infinity or Chebyshev) norm, under the
	c	 assumption that the boundaries are periodic or zero.  Add the
	c	 boundaries in with half weight (quarter weight on the edges
	c	 and eighth weight at the corners) for inhomogeneous boundaries.
	c-------------------------------------------------------------------*/

	double s, a;
	double tmp;
  int i3, i2, i1, n;
	dash::Array<double> p_s(dash::size());
	dash::Array<double> p_a(dash::size());

	auto pattern = r.pattern();
	auto local_beg_gidx = pattern.coords(pattern.global(0));
  auto local_end_gidx = pattern.coords(pattern.global(pattern.local_size()-1));

	int start = 0;
	if(local_beg_gidx[0] == 0) start++;

	int end = r.local.extent(0);
	if(local_end_gidx[0] == n3-1) end--;

	n = nx*ny*nz;

	for (i3 = start; i3 < end; i3++) {
    	for (i2 = 1; i2 < n2-1; i2++) {
            for (i1 = 1; i1 < n1-1; i1++) {
        		p_s.local[0] = p_s.local[0] + r.local(i3,i2,i1) * r.local(i3,i2,i1);
        		tmp = fabs(r.local(i3,i2,i1));
        		if (tmp > p_a.local[0]) p_a.local[0] = tmp;
        		}
    	}
  }

	dash::barrier();
	s = dash::reduce(p_s.begin(), p_s.end(), 0.0, std::plus<double>());
	a = (double) (*dash::max_element(p_a.begin(), p_a.end()));
	if (a > *rnmu) *rnmu = a;

	*rnm2 = sqrt(s/(double)n);

}

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void rep_nrm( dash::NArray<double, 3> &u, int n1, int n2, int n3, char *title, int kk) {
	/*--------------------------------------------------------------------
	c-------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
	c	 report on norm
	c-------------------------------------------------------------------*/

	double rnm2, rnmu;
	norm2u3(u,n1,n2,n3,&rnm2,&rnmu,nx[kk],ny[kk],nz[kk]);
	if(dash::myid() == 0) printf(" Level%2d in %8s: norms =%21.14e%21.14e\n", kk, title, rnm2, rnmu);
}

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void comm3( dash::NArray<double, 3> &u, int n1, int n2, int n3) {
	/*--------------------------------------------------------------------
	c-------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
	c	 comm3 organizes the communication on all borders
	c-------------------------------------------------------------------*/
	auto pattern = u.pattern();
	auto local_beg_gidx = pattern.coords(pattern.global(0));
  auto local_end_gidx = pattern.coords(pattern.global(pattern.local_size()-1));

	int start = 0;
	if(local_beg_gidx[0] == 0) start++;

	int end = u.local.extent(0);
	if(local_end_gidx[0] == n3-1) end--;

	// axis = 1

	for(int i3 = start; i3 < end; i3++) {
		for (int i2 = 1; i2 < n2-1; i2++) {
			u.local(i3,i2,n1-1) = u.local(i3,i2,1);
			u.local(i3,i2,0) = u.local(i3,i2,n1-2);
		}
	}

	// axis = 2

	for(int i3 = start; i3 < end; i3++) {
		for (int i1 = 0; i1 < n1; i1++) {
			u.local(i3,n2-1,i1) = u.local(i3,1,i1);
			u.local(i3,0,i1) = u.local(i3,n2-2,i1);
		}
	}

	dash::barrier();
	// axis = 3
	if(u(u.extent(0)-1,0,0).is_local()) {
		dash::copy(u.begin()+n2*n1*1, u.begin()+n2*n1*2, u.lbegin()+n2*n1*end);
	}

	if(u(0,0,0).is_local()) {
		dash::copy(u.begin()+n2*n1*(n3-2), u.begin()+n2*n1*(n3-1), u.lbegin());
	}

	dash::barrier();
}

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void zran3( dash::NArray<double, 3> &z, int n1, int n2, int n3, int nx, int ny, int k) {

	/*--------------------------------------------------------------------
	c-------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
	c	 zran3  loads +1 at ten randomly chosen points,
	c	 loads -1 at a different ten random points,
	c	 and zero elsewhere.
	c-------------------------------------------------------------------*/

	#define MM	10
	#define	A	pow(5.0,13)
	#define	X	314159265.e0

	int i0, m0, m1;
	/*int i1, i2, i3, d1, e1, e2, e3;*/
	int i1, i2, i3, d1, e2, e3;
	double xx, x0, x1, a1, a2, ai;

	double ten[MM][2], best;
	int i, j1[MM][2], j2[MM][2], j3[MM][2];


	/*double rdummy;*/

	a1 = power( A, nx );
	a2 = power( A, nx*ny );

	zero3(z,n1,n2,n3);

	i = is1-1+nx*(is2-1+ny*(is3-1));

	ai = power( A, i );
	d1 = ie1 - is1 + 1;
	/*e1 = ie1 - is1 + 2;*/
	e2 = ie2 - is2 + 2;
	e3 = ie3 - is3 + 2;
	x0 = X;
	/*rdummy = */randlc( &x0, ai );
	auto pattern = z.pattern();
	auto local_beg_gidx = pattern.coords(pattern.global(0));
	auto local_end_gidx = pattern.coords(pattern.global(pattern.local_size()-1));

	int start = 0;
	if(local_beg_gidx[0] == 0) start++;

	int end = z.local.extent(0);
	if(local_end_gidx[0] == z.extent(0)-1) end--;

	for(int j = 0; j < (int) local_beg_gidx[0]-1; j++){
		randlc( &x0, a2 );
	}

	for (i3 = start; i3 < end; i3++) {
		x1 = x0;
		for (i2 = 1; i2 < e2; i2++) {
			xx = x1;
			double tmp[d1+1];
			vranlc( d1, &xx, A, tmp);
			for(int i = 1; i <= d1; i++) z.local(i3,i2,i) = tmp[i];
			// vranlc( d1, &xx, A, (double *) &(z[i3][i2][0]));
			/*rdummy = */randlc( &x1, a1 );
		}
		/*rdummy = */randlc( &x0, a2 );
	}

	dash::barrier();
	if(dash::myid() == 0) {

	/*--------------------------------------------------------------------
	c	 call comm3(z,n1,n2,n3)
	c	 call showall(z,n1,n2,n3)
	c-------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
	c	 each processor looks for twenty candidates
	c-------------------------------------------------------------------*/

	for (i = 0; i < MM; i++) {
		ten[i][1] = 0.0;
		j1[i][1] = 0;
		j2[i][1] = 0;
		j3[i][1] = 0;
		ten[i][0] = 1.0;
		j1[i][0] = 0;
		j2[i][0] = 0;
		j3[i][0] = 0;
	}
	double current;
	for (i3 = 1; i3 < n3-1; i3++) {
		for (i2 = 1; i2 < n2-1; i2++) {
			for (i1 = 1; i1 < n1-1; i1++) {
				current = z(i3,i2,i1);
				if ( current > ten[0][1] ) {
					ten[0][1] = current;
					j1[0][1] = i1;
					j2[0][1] = i2;
					j3[0][1] = i3;
					bubble( ten, j1, j2, j3, MM, 1 );
				}
				if ( current < ten[0][0] ) {
					ten[0][0] = current;
					j1[0][0] = i1;
					j2[0][0] = i2;
					j3[0][0] = i3;
					bubble( ten, j1, j2, j3, MM, 0 );
				}
			}
		}
	}
	/*--------------------------------------------------------------------
	c	 Now which of these are globally best?
	c-------------------------------------------------------------------*/
	// i1 = MM - 1;
	// i0 = MM - 1;
	// int jg[4][MM][2];
	// for (i = MM - 1 ; i >= 0; i--) {
	// 	best = z[j3[i1][1]][j2[i1][1]][j1[i1][1]];
	// 	if (best == z[j3[i1][1]][j2[i1][1]][j1[i1][1]]) {
	// 		jg[0][i][1] = 0;
	// 		jg[1][i][1] = is1 - 1 + j1[i1][1];
	// 		jg[2][i][1] = is2 - 1 + j2[i1][1];
	// 		jg[3][i][1] = is3 - 1 + j3[i1][1];
	// 		i1 = i1-1;
	// 	} else {
	// 		jg[0][i][1] = 0;
	// 		jg[1][i][1] = 0;
	// 		jg[2][i][1] = 0;
	// 		jg[3][i][1] = 0;
	// 	}
	// 	ten[i][1] = best;
	// 	best = z[j3[i0][0]][j2[i0][0]][j1[i0][0]];
	// 	if (best == z[j3[i0][0]][j2[i0][0]][j1[i0][0]]) {
	// 		jg[0][i][0] = 0;
	// 		jg[1][i][0] = is1 - 1 + j1[i0][0];
	// 		jg[2][i][0] = is2 - 1 + j2[i0][0];
	// 		jg[3][i][0] = is3 - 1 + j3[i0][0];
	// 		i0 = i0-1;
	// 	} else {
	// 		jg[0][i][0] = 0;
	// 		jg[1][i][0] = 0;
	// 		jg[2][i][0] = 0;
	// 		jg[3][i][0] = 0;
	// 	}
	// 	ten[i][0] = best;
	// }
	// m1 = i1+1;
	// m0 = i0+1;

	/* printf(" negative charges at");
	for (i = 0; i < MM; i++) {
		if (i%5 == 0) printf("\n");
		printf(" (%3d,%3d,%3d)", jg[1][i][0], jg[2][i][0], jg[3][i][0]);
	}
	printf("\n positive charges at");
	for (i = 0; i < MM; i++) {
		if (i%5 == 0) printf("\n");
		printf(" (%3d,%3d,%3d)", jg[1][i][1], jg[2][i][1], jg[3][i][1]);
	}
	printf("\n small random numbers were\n");
	for (i = MM-1; i >= 0; i--) {
		printf(" %15.8e", ten[i][0]);
	}
	printf("\n and they were found on processor number\n");
	for (i = MM-1; i >= 0; i--) {
		printf(" %4d", jg[0][i][0]);
	}
	printf("\n large random numbers were\n");
	for (i = MM-1; i >= 0; i--) {
		printf(" %15.8e", ten[i][1]);
	}
	printf("\n and they were found on processor number\n");
	for (i = MM-1; i >= 0; i--) {
		printf(" %4d", jg[0][i][1]);
	}
	printf("\n");*/
}
	dash::barrier();
	zero3(z, n1, n2, n3);
	dash::barrier();
	if(dash::myid() == 0) {

		for (i = MM-1; i >= m0; i--) {
			z[j3[i][0]][j2[i][0]][j1[i][0]] = -1.0;
		}
		for (i = MM-1; i >= m1; i--) {
			z[j3[i][1]][j2[i][1]][j1[i][1]] = 1.0;
		}
	}
	dash::barrier();
	comm3(z,n1,n2,n3);

	/*--------------------------------------------------------------------
	c	 call showall(z,n1,n2,n3)
	c-------------------------------------------------------------------*/

}

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void showall( dash::NArray<double, 3> &z, int n1, int n2, int n3) {
	if(dash::myid() == 0) {
		/*--------------------------------------------------------------------
		c-------------------------------------------------------------------*/

		int i1,i2,i3;
		int m1, m2, m3;

		m1 = min(n1,18);
		m2 = min(n2,14);
		m3 = min(n3,18);

		printf("\n");
		for (i3 = 0; i3 < m3; i3++) {
			for (i1 = 0; i1 < m1; i1++) {
				for (i2 = 0; i2 < m2; i2++) {
				  printf("%6.3f", (double) z[i3][i2][i1]);
				}
				printf("\n");
			}
			printf(" - - - - - - - \n");
		}
		printf("\n");
	}
}

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static double power( double a, int n ) {

	/*--------------------------------------------------------------------
	c-------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
	c	 power  raises an integer, disguised as a double
	c	 precision real, to an integer power
	c-------------------------------------------------------------------*/
	double aj;
	int nj;
	/* double rdummy;*/
	double power;

	power = 1.0;
	nj = n;
	aj = a;

	while (nj != 0) {
		if( (nj%2) == 1 ) /*rdummy =  */randlc( &power, aj );
		/*rdummy = */randlc( &aj, aj );
		nj = nj/2;
	}

	return (power);
}

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void bubble( double ten[M][2], int j1[M][2], int j2[M][2], int j3[M][2], int m, int ind ) {

	/*--------------------------------------------------------------------
	c-------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
	c	 bubble		does a bubble sort in direction dir
	c-------------------------------------------------------------------*/
	double temp;
	int i, j_temp;
	if ( ind == 1 ) {
		for (i = 0; i < m-1; i++) {
			if ( ten[i][ind] > ten[i+1][ind] ) {
				temp = ten[i+1][ind];
				ten[i+1][ind] = ten[i][ind];
				ten[i][ind] = temp;

				j_temp = j1[i+1][ind];
				j1[i+1][ind] = j1[i][ind];
				j1[i][ind] = j_temp;

				j_temp = j2[i+1][ind];
				j2[i+1][ind] = j2[i][ind];
				j2[i][ind] = j_temp;

				j_temp = j3[i+1][ind];
				j3[i+1][ind] = j3[i][ind];
				j3[i][ind] = j_temp;
			} else {
				return;
			}
		}
		} else {
		for (i = 0; i < m-1; i++) {
			if ( ten[i][ind] < ten[i+1][ind]){

				temp = ten[i+1][ind];
				ten[i+1][ind] = ten[i][ind];
				ten[i][ind] = temp;

				j_temp = j1[i+1][ind];
				j1[i+1][ind] = j1[i][ind];
				j1[i][ind] = j_temp;

				j_temp = j2[i+1][ind];
				j2[i+1][ind] = j2[i][ind];
				j2[i][ind] = j_temp;

				j_temp = j3[i+1][ind];
				j3[i+1][ind] = j3[i][ind];
				j3[i][ind] = j_temp;
			} else {
				return;
			}
		}
	}
}

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

static void zero3( dash::NArray<double, 3> &z, int n1, int n2, int n3) {
	/*--------------------------------------------------------------------
	c-------------------------------------------------------------------*/

	// for(int i3 = 0; i3 < z.local.extent(0); i3++){
	// 	for (int i2 = 0; i2 < n2; i2++) {
	// 		for (int i1 = 0; i1 < n1; i1++) {
	// 			z.local(i3,i2,i1) = 0.0;
	// 		}
	// 	}
	// }
	dash::fill(z.begin(), z.end(), 0.0);
}

/*---- end of program ------------------------------------------------*/
