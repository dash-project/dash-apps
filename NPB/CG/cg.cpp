/*--------------------------------------------------------------------

	Information on NAS Parallel Benchmarks is available at:

	http://www.nas.nasa.gov/Software/NPB/

	Authors: M. Yarrow
	C. Kuszmaul

	STL version:
	Nicco Mietzsch <nicco.mietzsch@campus.lmu.de>

	CPP and OpenMP version:
			Dalvan Griebler <dalvangriebler@gmail.com>
			Júnior Löff <loffjh@gmail.com>

--------------------------------------------------------------------*/

/*
c---------------------------------------------------------------------
c  Note: please observe that in the routine conj_grad three
c  implementations of the sparse matrix-vector multiply have
c  been supplied.  The default matrix-vector multiply is not
c  loop unrolled.  The alternate implementations are unrolled
c  to a depth of 2 and unrolled to a depth of 8.  Please
c  experiment with these to find the fastest for your particular
c  architecture.  If reporting timing results, any of these three may
c  be used without penalty.
c---------------------------------------------------------------------
*/
#include <libdash.h>

#include <iostream>
#include "npbparams.hpp"
#include "npb-CPP.hpp"

#define	NZ	NA*(NONZER+1)*(NONZER+1)+NA*(NONZER+2)

#define TIMER_ENABLED FALSE

/* global variables */

/* common /partit_size/ */
static int naa;
static int nzz;
static int firstrow;
static int lastrow;
static int firstcol;
static int lastcol;

/* common /urando/ */
static double amult;
static double tran;

/* function declarations */
static void conj_grad (dash::Array<int, int, dash::CSRPattern<1>> &colidx, dash::Array<int> &rowstr, dash::Array<double> &x,
	dash::Array<double> &z, dash::Array<double, int, dash::CSRPattern<1>> &al, dash::Array<double> &p, dash::Array<double> &q,
	dash::Array<double> &r, dash::Array<double> &w, double *rnorm);
static void makea(int n, int nz, double a[], int colidx[], int rowstr[],
	int nonzer, int firstrow, int lastrow, int firstcol,
	int lastcol, double rcond, int arow[], int acol[],
	double aelt[], double v[], int iv[], double shift );
static void sparse(double a[], int colidx[], int rowstr[], int n,
	int arow[], int acol[], double aelt[],
	int firstrow, int lastrow,
	double x[], boolean mark[], int nzloc[], int nnza);
static void sprnvc(int n, int nz, double v[], int iv[], int nzloc[], int mark[]);
static int icnvrt(double x, int ipwr2);
static void vecset(int n, double v[], int iv[], int *nzv, int i, double val);

/*--------------------------------------------------------------------
      program cg
--------------------------------------------------------------------*/

int main(int argc, char **argv)
{

	dash::init(&argc, &argv);

	dash::Array<double> a_global(NZ+1);
	dash::Array<int> rowstr_global(NA+1+1);
	dash::Array<int> colidx_global(NZ+1);

	dash::Array<double> x(NA+1);	/* x[1:NA] */
	dash::Array<double> z(NA+1);	/* z[1:NA] */
	dash::Array<double> p(NA+1);	/* p[1:NA] */
	dash::Array<double> q(NA+1);	/* q[1:NA] */
	dash::Array<double> r(NA+1);	/* r[1:NA] */
	dash::Array<double> w(NA+1);	/* w[1:NA] */


	int i, j, k, it;
	int nthreads = dash::size();
	double zeta;
	double rnorm;
	double norm_temp11;
	double norm_temp12;
	double t, mflops;
	char class_npb;
	boolean verified;
	double zeta_verify_value, epsilon;

	firstrow = 1;
	lastrow  = NA;
	firstcol = 1;
	lastcol  = NA;

	if (NA == 1400 && NONZER == 7 && NITER == 15 && SHIFT == 10.0) {
		class_npb = 'S';
		zeta_verify_value = 8.5971775078648;
	} else if (NA == 7000 && NONZER == 8 && NITER == 15 && SHIFT == 12.0) {
		class_npb = 'W';
		zeta_verify_value = 10.362595087124;
	} else if (NA == 14000 && NONZER == 11 && NITER == 15 && SHIFT == 20.0) {
		class_npb = 'A';
		zeta_verify_value = 17.130235054029;
	} else if (NA == 75000 && NONZER == 13 && NITER == 75 && SHIFT == 60.0) {
		class_npb = 'B';
		zeta_verify_value = 22.712745482631;
	} else if (NA == 150000 && NONZER == 15 && NITER == 75 && SHIFT == 110.0) {
		class_npb = 'C';
		zeta_verify_value = 28.973605592845;
	} else {
		class_npb = 'U';
	}

	if (0 == dash::myid()) {

		printf("\n\n NAS Parallel Benchmarks 4.0 C++ DASH version"" - CG Benchmark\n");
		printf("\n\n Developed by: Dalvan Griebler <dalvan.griebler@acad.pucrs.br>\n");
		printf("\n\n DASH version by: Nicco Mietzsch <nicco.mietzsch@campus.lmu.de>\n");
		printf(" Size: %10d\n", NA);
		printf(" Iterations: %5d\n", NITER);


		naa = NA;
		nzz = NZ;

		/*--------------------------------------------------------------------
		c  Initialize random number generator
		c-------------------------------------------------------------------*/
		tran    = 314159265.0;
		amult   = 1220703125.0;
		zeta    = randlc( &tran, amult );

		/*--------------------------------------------------------------------
		c
		c-------------------------------------------------------------------*/
		static int colidx_init[NZ+1];	/* colidx[1:NZ] */
		static int rowstr_init[NA+1+1];	/* rowstr[1:NA+1] */

		static double a_init[NZ+1];		/* a[1:NZ] */



		static int iv[2*NA+1+1];	/* iv[1:2*NA+1] */
		static int arow[NZ+1];		/* arow[1:NZ] */
		static int acol[NZ+1];		/* acol[1:NZ] */

		static double v[NA+1+1];	/* v[1:NA+1] */
		static double aelt[NZ+1];	/* aelt[1:NZ] */


		makea(naa, nzz, a_init, colidx_init, rowstr_init, NONZER,
		firstrow, lastrow, firstcol, lastcol,
		RCOND, arow, acol, aelt, v, iv, SHIFT);

		for(i = 0; i < NZ+1; ++i) colidx_global[i] = colidx_init[i];
		for(i = 0; i < NA+1+1; ++i) rowstr_global[i] = rowstr_init[i];
		for(i = 0; i < NZ+1; ++i) a_global[i] = a_init[i];


		/*--------------------------------------------------------------------
		c  set starting vector to (1, 1, .... 1)
		c-------------------------------------------------------------------*/

		for (i = 1; i <= NA; i++) {
			x[i] = 1.0;
		}

	}
	dash::barrier();

	//setup local arrays
	// first, setup rowstr

	int elem_per_unit = ceil(((double) NA+1)/dash::size());

	dash::Array<int> rowstr((elem_per_unit+1)*dash::size());

	for(int i = 0; i < elem_per_unit+1; i++) {
			rowstr.local[i] = rowstr_global[min(dash::myid()*elem_per_unit + i, NA+1)];
	}


	//now setup a and colidx

	typedef dash::CSRPattern<1>           pattern_t;
  typedef int                           index_t;
	typedef typename pattern_t::size_type extent_t;

	std::vector<extent_t> local_sizes;
	dash::Array<int> my_elem_count(dash::size());

	my_elem_count.local[0] = rowstr.local[elem_per_unit] - rowstr.local[0];

	for (int unit_idx = 0; unit_idx < dash::size(); ++unit_idx) {
    local_sizes.push_back(my_elem_count[unit_idx]);
  }

	pattern_t pattern(local_sizes);

	dash::Array<double, index_t, pattern_t> a(pattern);
	dash::Array<int, index_t, pattern_t> colidx(pattern);

	int offset = 0;
	for(int k = 0; k < dash::myid(); ++k) offset += my_elem_count[k];

	for(int i = 0; i < my_elem_count.local[0]; ++i) {
		a.local[i] = a_global[offset + i];
		colidx.local[i] = colidx_global[offset + i];
	}

	for(int i = 0; i < elem_per_unit+1; ++i) {
		rowstr.local[i] -= offset;
	}

	dash::barrier();

	zeta = 0.0;

	/*-------------------------------------------------------------------
	c---->
	c  Do one iteration untimed to init all code and data page tables
	c---->                    (then reinit, start timing, to niter its)
	c-------------------------------------------------------------------*/

	for (it = 1; it <= 1; it++) {

		/*--------------------------------------------------------------------
		c  The call to the conjugate gradient routine:
		c-------------------------------------------------------------------*/
		conj_grad (colidx, rowstr, x, z, a, p, q, r, w, &rnorm);

		/*--------------------------------------------------------------------
		c  zeta = shift + 1/(x.z)
		c  So, first: (x.z)
		c  Also, find norm of z
		c  So, first: (z.z)
		c-------------------------------------------------------------------*/
		dash::Array<double> norm_temp11(dash::size());

		auto it2 = z.lbegin();
		for( auto it1 = x.lbegin(); it1 != x.lend(); ++it1) {
			norm_temp11.local[0] += *it1 * *it2;
			it2++;
		}

		norm_temp11.local[0] = dash::reduce(norm_temp11.begin(), norm_temp11.end(), 0.0, std::plus<double>());

		dash::Array<double> norm_temp12(dash::size());

		for( auto it1 = z.lbegin(); it1 != z.lend(); ++it1) {
			norm_temp12.local[0] += *it1 * *it1;
		}

		norm_temp12.local[0] = dash::reduce(norm_temp12.begin(), norm_temp12.end(), 0.0, std::plus<double>());

		norm_temp12.local[0] = 1.0 / sqrt( norm_temp12.local[0] );

		/*--------------------------------------------------------------------
		c  Normalize z to obtain x
		c-------------------------------------------------------------------*/

		it2 = z.lbegin();
		for( auto it1 = x.lbegin(); it1 != x.lend(); ++it1) {
			*it1 = *it2 * norm_temp12.local[0];
			it2++;
		}

	} /* end of do one iteration untimed */

	/*--------------------------------------------------------------------
	c  set starting vector to (1, 1, .... 1)
	c-------------------------------------------------------------------*/
	if (0 == dash::myid()) {

		for (i = 1; i <= NA; i++) {
			x[i] = 1.0;
		}

	}

	zeta = 0.0;

	dash::barrier();

	if (0 == dash::myid()) {

		timer_clear( 1 );
		timer_clear( 2 );

		timer_start( 1 );

	}

	/*--------------------------------------------------------------------
	c---->
	c  Main Iteration for inverse power method
	c---->
	c-------------------------------------------------------------------*/


	for (it = 1; it <= NITER; it++) {

		/*--------------------------------------------------------------------
		c  The call to the conjugate gradient routine:
		c-------------------------------------------------------------------*/
		conj_grad(colidx, rowstr, x, z, a, p, q, r, w, &rnorm);
		/*--------------------------------------------------------------------
		c  zeta = shift + 1/(x.z)
		c  So, first: (x.z)
		c  Also, find norm of z
		c  So, first: (z.z)
		c-------------------------------------------------------------------*/

		if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_start(2);

		dash::Array<double> norm_temp11(dash::size());

		auto it2 = z.lbegin();
		for( auto it1 = x.lbegin(); it1 != x.lend(); ++it1) {
			norm_temp11.local[0] += *it1 * *it2;
			it2++;
		}

		norm_temp11.local[0] = dash::reduce(norm_temp11.begin(), norm_temp11.end(), 0.0, std::plus<double>());

		dash::Array<double> norm_temp12(dash::size());

		for( auto it1 = z.lbegin(); it1 != z.lend(); ++it1) {
			norm_temp12.local[0] += *it1 * *it1;
		}

		norm_temp12.local[0] = dash::reduce(norm_temp12.begin(), norm_temp12.end(), 0.0, std::plus<double>());

		if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_stop(2);

		norm_temp12.local[0] = 1.0 / sqrt( norm_temp12.local[0] );

		if(0 == dash::myid()) {

			zeta = SHIFT + 1.0 / norm_temp11.local[0];

			if( it == 1 ) {
				printf("   iteration           ||r||                 zeta\n");
			}

			printf("    %5d       %20.14e%20.13e\n", it, rnorm, zeta);

		}

		/*--------------------------------------------------------------------
		c  Normalize z to obtain x
		c-------------------------------------------------------------------*/

		if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_start(2);

		it2 = z.lbegin();
		for( auto it1 = x.lbegin(); it1 != x.lend(); ++it1) {
			*it1 = *it2 * norm_temp12.local[0];
			it2++;
		}

		if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_stop(2);

		/* end of main iter inv pow meth */

	}

	dash::barrier();

	if (0 == dash::myid()) {


		timer_stop( 1 );
		/*--------------------------------------------------------------------
		c  End of timed section
		c-------------------------------------------------------------------*/

		t = timer_read( 1 );

		printf(" Benchmark completed\n");

		epsilon = 1.0e-10;
		if (class_npb != 'U') {
			if (fabs(zeta - zeta_verify_value) <= epsilon) {
				verified = TRUE;
				printf(" VERIFICATION SUCCESSFUL\n");
				printf(" Zeta is    %20.12e\n", zeta);
				printf(" Error is   %20.12e\n", zeta - zeta_verify_value);
			} else {
				verified = FALSE;
				printf(" VERIFICATION FAILED\n");
				printf(" Zeta                %20.12e\n", zeta);
				printf(" The correct zeta is %20.12e\n", zeta_verify_value);
			}
		} else {
			verified = FALSE;
			printf(" Problem size unknown\n");
			printf(" NO VERIFICATION PERFORMED\n");
		}
		if ( t != 0.0 ) {
			mflops = (2.0*NITER*NA)
			* (3.0+(NONZER*(NONZER+1)) + 25.0*(5.0+(NONZER*(NONZER+1))) + 3.0 )
			/ t / 1000000.0;
		} else {
			mflops = 0.0;
		}
		c_print_results((char*)"CG", class_npb, NA, 0, 0, NITER, nthreads, t, mflops, (char*)"          floating point", verified, (char*)NPBVERSION, (char*)COMPILETIME, (char*)CS1, (char*)CS2, (char*)CS3, (char*)CS4, (char*)CS5, (char*)CS6, (char*)CS7);

	}

	return 0;
}

/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/
static void conj_grad (
	dash::Array<int, int, dash::CSRPattern<1>> &colidx,	/* colidx[1:nzz] */
	dash::Array<int> &rowstr,	/* rowstr[1:naa+1] */
	dash::Array<double> &x,	/* x[*] */
	dash::Array<double> &z,	/* z[*] */
	dash::Array<double, int, dash::CSRPattern<1>> &a,	/* a[1:nzz] */
	dash::Array<double> &p,	/* p[*] */
	dash::Array<double> &q,	/* q[*] */
	dash::Array<double> &r,	/* r[*] */
	dash::Array<double> &w,	/* w[*] */
	double *rnorm )
/*--------------------------------------------------------------------
c-------------------------------------------------------------------*/

/*---------------------------------------------------------------------
c  Floaging point arrays here are named as in NPB1 spec discussion of
c  CG algorithm
c---------------------------------------------------------------------*/
{
	//static double d, sum, rho, rho0, alpha, beta;
	//double d, sum, rho, rho0, alpha, beta;
	dash::Array<double> d(dash::size());
	dash::Array<double> sum(dash::size());
	dash::Array<double> rho(dash::size());
	dash::Array<double> rho0(dash::size());
	dash::Array<double> alpha(dash::size());
	dash::Array<double> beta(dash::size());
	int j;
	int cgit, cgitmax = 25;


	rho.local[0] = 0.0;

	int elem_per_unit = ceil(((double) NA+1)/dash::size());

	/*--------------------------------------------------------------------
	c  Initialize the CG algorithm:
	c-------------------------------------------------------------------*/

	if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_start(2);

	dash::fill(q.begin(), q.end(), 0.0);
	dash::fill(z.begin(), z.end(), 0.0);
	dash::fill(w.begin(), w.end(), 0.0);

	std::copy(x.lbegin(), x.lend(), r.lbegin());
	std::copy(r.lbegin(), r.lend(), p.lbegin());

	if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_stop(2);

	/*--------------------------------------------------------------------
	c  rho = r.r
	c  Now, obtain the norm of r: First, sum squares of r elements locally...
	c-------------------------------------------------------------------*/

	if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_start(2);

	for( auto it1 = r.lbegin(); it1 != r.lend(); it1++) {
		rho.local[0] = rho.local[0] +  *it1 * *it1;
	}

	rho.local[0] = dash::reduce(rho.begin(), rho.end(), 0.0, std::plus<double>());

	if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_stop(2);

	/*--------------------------------------------------------------------
	c---->
	c  The conj grad iteration loop
	c---->
	c-------------------------------------------------------------------*/

	for (cgit = 1; cgit <= cgitmax; cgit++) {

		rho0.local[0] = rho.local[0];
		d.local[0] = 0.0;
		rho.local[0] = 0.0;


		/*--------------------------------------------------------------------
		c  q = A.p
		c  The partition submatrix-vector multiply: use workspace w
		c---------------------------------------------------------------------
		C
		C  NOTE: this version of the multiply is actually (slightly: maybe %5)
		C        faster on the sp2 on 16 nodes than is the unrolled-by-2 version
		C        below.   On the Cray t3d, the reverse is true, i.e., the
		C        unrolled-by-two version is some 10% faster.
		C        The unrolled-by-8 version below is significantly faster
		C        on the Cray t3d - overall speed of code is 1.5 times faster.
		*/

		double p_local[p.size()];
		dash::copy(p.begin(), p.end(), &p_local[0]);

		if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_start(2);

		// rolled version

		for(int j = 0; j < elem_per_unit; j++) {
			double sum = 0.0;
			for (int k = rowstr.local[j]; k < rowstr.local[j+1]; k++) {
				sum = sum + a.local[k]*p_local[colidx.local[k]];
			}
			w.local[j] = sum;
		}


		//unrolled-by-two version
		/*
		for(int j = 0; j < elem_per_unit; j++) {
			int iresidue;
			double sum1, sum2;
			int i = rowstr.local[j];
			iresidue = (rowstr.local[j+1]-i) % 2;
			sum1 = 0.0;
			sum2 = 0.0;
			if (iresidue == 1) sum1 = sum1 + a.local[i]*p_local[colidx.local[i]];
			for (int k = i+iresidue; k <= rowstr.local[j+1]-2; k += 2) {
				sum1 = sum1 + a.local[k]   * p_local[colidx.local[k]];
				sum2 = sum2 + a.local[k+1] * p_local[colidx.local[k+1]];
			}
			w.local[j] = sum1 + sum2;
		}*/


		//unrolled-by-8 version
		/*
		for(int j = 0; j < elem_per_unit; j++) {
			int iresidue, k;
			int i = rowstr.local[j];
			iresidue = (rowstr.local[j+1]-i) % 8;
			double sum = 0.0;
			for (k = i; k <= i+iresidue-1; k++) {
				sum = sum +  a.local[k] * p_local[colidx.local[k]];
			}
			for (k = i+iresidue; k <= rowstr.local[j+1]-8; k += 8) {
				sum = sum
				+ a.local[k  ] * p_local[colidx.local[k  ]]
				+ a.local[k+1] * p_local[colidx.local[k+1]]
				+ a.local[k+2] * p_local[colidx.local[k+2]]
				+ a.local[k+3] * p_local[colidx.local[k+3]]
				+ a.local[k+4] * p_local[colidx.local[k+4]]
				+ a.local[k+5] * p_local[colidx.local[k+5]]
				+ a.local[k+6] * p_local[colidx.local[k+6]]
				+ a.local[k+7] * p_local[colidx.local[k+7]];
			}
			w.local[j] = sum;
		}*/


		if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_stop(2);

		if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_start(2);
		std::copy(w.lbegin(), w.lend(), q.lbegin());
		if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_stop(2);
		/*--------------------------------------------------------------------
		c  Clear w for reuse... Do we really need that? Commented it out for the moment.
		c-------------------------------------------------------------------*/

		if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_start(2);
		//dash::fill(w.begin(), w.end(), 0.0);
		if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_stop(2);
		/*--------------------------------------------------------------------
		c  Obtain p.q
		c-------------------------------------------------------------------*/

		if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_start(2);

		auto it2 = p.lbegin();
		for( auto it1 = q.lbegin(); it1 != q.lend(); ++it1) {
			d.local[0] += *it1 * *it2;
			it2++;
		}

		d.local[0] = dash::reduce(d.begin(), d.end(), 0.0, std::plus<double>());

		if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_stop(2);
		/*--------------------------------------------------------------------
		c  Obtain alpha = rho / (p.q)
		c-------------------------------------------------------------------*/

		alpha.local[0] = rho0.local[0] / d.local[0];

		/*--------------------------------------------------------------------
		c  Save a temporary of rho
		c-------------------------------------------------------------------*/
			/*	rho0 = rho;*/

		/*---------------------------------------------------------------------
		c  Obtain z = z + alpha*p
		c  and    r = r - alpha*q
		c---------------------------------------------------------------------*/

		if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_start(2);

		it2 = p.lbegin();
		for( auto it1 = z.lbegin(); it1 != z.lend(); ++it1) {
			*it1 = *it1 + alpha.local[0]* *it2;
			it2++;
		}

		it2 = q.lbegin();
		for( auto it1 = r.lbegin(); it1 != r.lend(); ++it1) {
			*it1 = *it1 - alpha.local[0]* *it2;
			it2++;
		}

		if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_stop(2);
		/*---------------------------------------------------------------------
		c  rho = r.r
		c  Now, obtain the norm of r: First, sum squares of r elements locally...
		c---------------------------------------------------------------------*/

		if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_start(2);

		for( auto it1 = r.lbegin(); it1 != r.lend(); ++it1) {
			rho.local[0] += *it1 * *it1;
		}

		rho.local[0] = dash::reduce(rho.begin(), rho.end(), 0.0, std::plus<double>());

		if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_stop(2);
		/*--------------------------------------------------------------------
		c  Obtain beta:
		c-------------------------------------------------------------------*/

		beta.local[0] = rho.local[0] / rho0.local[0];

		/*--------------------------------------------------------------------
		c  p = r + beta*p
		c-------------------------------------------------------------------*/

		if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_start(2);

		it2 = r.lbegin();
		for( auto it1 = p.lbegin(); it1 != p.lend(); ++it1) {
			*it1 = *it2 + beta.local[0]* *it1;
			it2++;
		}
		dash::barrier();

		if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_stop(2);
	} /* end of do cgit=1,cgitmax */

	/*---------------------------------------------------------------------
	c  Compute residual norm explicitly:  ||r|| = ||x - A.z||
	c  First, form A.z
	c  The partition submatrix-vector multiply
	c---------------------------------------------------------------------*/

	sum.local[0] = 0.0;
	if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_start(2);

	double z_local[p.size()];
	dash::copy(z.begin(), z.end(), &z_local[0]);

	for(int j = 0; j < elem_per_unit; j++) {
		double sum = 0.0;
		for (int k = rowstr.local[j]; k < rowstr.local[j+1]; k++) {
			sum = sum + a.local[k]*z_local[colidx.local[k]];
		}
		w.local[j] = sum;
	}

	if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_stop(2);

	if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_start(2);
	std::copy(w.lbegin(), w.lend(), r.lbegin());
	if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_stop(2);
	/*--------------------------------------------------------------------
	c  At this point, r contains A.z
	c-------------------------------------------------------------------*/
	if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_start(2);

	auto it2 = r.lbegin();
	for( auto it1 = x.lbegin(); it1 != x.lend(); ++it1) {
		sum.local[0] += (*it1 - *it2)*(*it1 - *it2);
		it2++;
	}

	sum.local[0] = dash::reduce(sum.begin(), sum.end(), 0.0, std::plus<double>());

	if(TIMER_ENABLED == TRUE && 0 == dash::myid()) timer_stop(2);

	(*rnorm) = sqrt(sum.local[0]);
}

/*---------------------------------------------------------------------
c       generate the test problem for benchmark 6
c       makea generates a sparse matrix with a
c       prescribed sparsity distribution
c
c       parameter    type        usage
c
c       input
c
c       n            i           number of cols/rows of matrix
c       nz           i           nonzeros as declared array size
c       rcond        r*8         condition number
c       shift        r*8         main diagonal shift
c
c       output
c
c       a            r*8         array for nonzeros
c       colidx       i           col indices
c       rowstr       i           row pointers
c
c       workspace
c
c       iv, arow, acol i
c       v, aelt        r*8
c---------------------------------------------------------------------*/
static void makea(
	int n,
	int nz,
	double a[],		/* a[1:nz] */
	int colidx[],	/* colidx[1:nz] */
	int rowstr[],	/* rowstr[1:n+1] */
	int nonzer,
	int firstrow,
	int lastrow,
	int firstcol,
	int lastcol,
	double rcond,
	int arow[],		/* arow[1:nz] */
	int acol[],		/* acol[1:nz] */
	double aelt[],	/* aelt[1:nz] */
	double v[],		/* v[1:n+1] */
	int iv[],		/* iv[1:2*n+1] */
	double shift )
{
	int i, nnza, iouter, ivelt, ivelt1, irow, nzv;

	/*--------------------------------------------------------------------
	c      nonzer is approximately  (int(sqrt(nnza /n)));
	c-------------------------------------------------------------------*/

	double size, ratio, scale;
	int jcol;

	size = 1.0;
	ratio = pow(rcond, (1.0 / (double)n));
	nnza = 0;

	/*---------------------------------------------------------------------
	c  Initialize colidx(n+1 .. 2n) to zero.
	c  Used by sprnvc to mark nonzero positions
	c---------------------------------------------------------------------*/

	for (i = 1; i <= n; i++) {
		colidx[n+i] = 0;
	}
	for (iouter = 1; iouter <= n; iouter++) {
		nzv = nonzer;
		sprnvc(n, nzv, v, iv, &(colidx[0]), &(colidx[n]));
		vecset(n, v, iv, &nzv, iouter, 0.5);
		for (ivelt = 1; ivelt <= nzv; ivelt++){
			jcol = iv[ivelt];
			if (jcol >= firstcol && jcol <= lastcol) {
				scale = size * v[ivelt];
				for (ivelt1 = 1; ivelt1 <= nzv; ivelt1++) {
					irow = iv[ivelt1];
					if (irow >= firstrow && irow <= lastrow) {
						nnza = nnza + 1;
						if (nnza > nz) {
							printf("Space for matrix elements exceeded in" " makea\n");
							printf("nnza, nzmax = %d, %d\n", nnza, nz);
							printf("iouter = %d\n", iouter);
							exit(1);
						}
						acol[nnza] = jcol;
						arow[nnza] = irow;
						aelt[nnza] = v[ivelt1] * scale;
					}
				}
			}
		}
		size = size * ratio;
	}

	/*---------------------------------------------------------------------
	c       ... add the identity * rcond to the generated matrix to bound
	c           the smallest eigenvalue from below by rcond
	c---------------------------------------------------------------------*/
	for (i = firstrow; i <= lastrow; i++) {
		if (i >= firstcol && i <= lastcol) {
			iouter = n + i;
			nnza = nnza + 1;
			if (nnza > nz) {
				printf("Space for matrix elements exceeded in makea\n");
				printf("nnza, nzmax = %d, %d\n", nnza, nz);
				printf("iouter = %d\n", iouter);
				exit(1);
			}
			acol[nnza] = i;
			arow[nnza] = i;
			aelt[nnza] = rcond - shift;
		}
	}

	/*---------------------------------------------------------------------
	c       ... make the sparse matrix from list of elements with duplicates
	c           (v and iv are used as  workspace)
	c---------------------------------------------------------------------*/
	sparse(a, colidx, rowstr, n, arow, acol, aelt, firstrow, lastrow, v, &(iv[0]), &(iv[n]), nnza);
}

/*---------------------------------------------------
c       generate a sparse matrix from a list of
c       [col, row, element] tri
c---------------------------------------------------*/
static void sparse(
	double a[],		/* a[1:*] */
	int colidx[],	/* colidx[1:*] */
	int rowstr[],	/* rowstr[1:*] */
	int n,
	int arow[],		/* arow[1:*] */
	int acol[],		/* acol[1:*] */
	double aelt[],	/* aelt[1:*] */
	int firstrow,
	int lastrow,
	double x[],		/* x[1:n] */
	boolean mark[],	/* mark[1:n] */
	int nzloc[],	/* nzloc[1:n] */
	int nnza)
/*---------------------------------------------------------------------
c       rows range from firstrow to lastrow
c       the rowstr pointers are defined for nrows = lastrow-firstrow+1 values
c---------------------------------------------------------------------*/
{
	int nrows;
	int i, j, jajp1, nza, k, nzrow;
	double xi;

	/*--------------------------------------------------------------------
	c    how many rows of result
	c-------------------------------------------------------------------*/
	nrows = lastrow - firstrow + 1;

	/*--------------------------------------------------------------------
	c     ...count the number of triples in each row
	c-------------------------------------------------------------------*/

	for (j = 1; j <= n; j++) {
		rowstr[j] = 0;
		mark[j] = FALSE;
	}
	rowstr[n+1] = 0;

	for (nza = 1; nza <= nnza; nza++) {
		j = (arow[nza] - firstrow + 1) + 1;
		rowstr[j] = rowstr[j] + 1;
	}
	rowstr[1] = 1;
	for (j = 2; j <= nrows+1; j++) {
		rowstr[j] = rowstr[j] + rowstr[j-1];
	}

	/*---------------------------------------------------------------------
	c     ... rowstr(j) now is the location of the first nonzero
	c           of row j of a
	c---------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
	c     ... do a bucket sort of the triples on the row index
	c-------------------------------------------------------------------*/
	for (nza = 1; nza <= nnza; nza++) {
		j = arow[nza] - firstrow + 1;
		k = rowstr[j];
		a[k] = aelt[nza];
		colidx[k] = acol[nza];
		rowstr[j] = rowstr[j] + 1;
	}

	/*--------------------------------------------------------------------
	c       ... rowstr(j) now points to the first element of row j+1
	c-------------------------------------------------------------------*/
	for (j = nrows; j >= 1; j--) {
		rowstr[j+1] = rowstr[j];
	}
	rowstr[1] = 1;

	/*--------------------------------------------------------------------
	c       ... generate the actual output rows by adding elements
	c-------------------------------------------------------------------*/
	nza = 0;

	for (i = 1; i <= n; i++) {
		x[i] = 0.0;
		mark[i] = FALSE;
	}

	jajp1 = rowstr[1];
	for (j = 1; j <= nrows; j++) {
		nzrow = 0;

		/*--------------------------------------------------------------------
		c          ...loop over the jth row of a
		c-------------------------------------------------------------------*/
		for (k = jajp1; k < rowstr[j+1]; k++) {
			i = colidx[k];
			x[i] = x[i] + a[k];
			if ( mark[i] == FALSE && x[i] != 0.0) {
				mark[i] = TRUE;
				nzrow = nzrow + 1;
				nzloc[nzrow] = i;
			}
		}

		/*--------------------------------------------------------------------
		c          ... extract the nonzeros of this row
		c-------------------------------------------------------------------*/
		for (k = 1; k <= nzrow; k++) {
			i = nzloc[k];
			mark[i] = FALSE;
			xi = x[i];
			x[i] = 0.0;
			if (xi != 0.0) {
				nza = nza + 1;
				a[nza] = xi;
				colidx[nza] = i;
			}
		}
		jajp1 = rowstr[j+1];
		rowstr[j+1] = nza + rowstr[1];
	}
}

/*---------------------------------------------------------------------
c       generate a sparse n-vector (v, iv)
c       having nzv nonzeros
c
c       mark(i) is set to 1 if position i is nonzero.
c       mark is all zero on entry and is reset to all zero before exit
c       this corrects a performance bug found by John G. Lewis, caused by
c       reinitialization of mark on every one of the n calls to sprnvc
---------------------------------------------------------------------*/
static void sprnvc(
	int n,
	int nz,
	double v[],		/* v[1:*] */
	int iv[],		/* iv[1:*] */
	int nzloc[],	/* nzloc[1:n] */
	int mark[] ) 	/* mark[1:n] */
{
	int nn1;
	int nzrow, nzv, ii, i;
	double vecelt, vecloc;

	nzv = 0;
	nzrow = 0;
	nn1 = 1;
	do {
		nn1 = 2 * nn1;
	} while (nn1 < n);

	/*--------------------------------------------------------------------
	c    nn1 is the smallest power of two not less than n
	c-------------------------------------------------------------------*/

	while (nzv < nz) {
		vecelt = randlc(&tran, amult);

		/*--------------------------------------------------------------------
		c   generate an integer between 1 and n in a portable manner
		c-------------------------------------------------------------------*/
		vecloc = randlc(&tran, amult);
		i = icnvrt(vecloc, nn1) + 1;
		if (i > n) continue;

		/*--------------------------------------------------------------------
		c  was this integer generated already?
		c-------------------------------------------------------------------*/
		if (mark[i] == 0) {
			mark[i] = 1;
			nzrow = nzrow + 1;
			nzloc[nzrow] = i;
			nzv = nzv + 1;
			v[nzv] = vecelt;
			iv[nzv] = i;
		}
	}

	for (ii = 1; ii <= nzrow; ii++) {
		i = nzloc[ii];
		mark[i] = 0;
	}
}

/*---------------------------------------------------------------------
* scale a double precision number x in (0,1) by a power of 2 and chop it
*---------------------------------------------------------------------*/
static int icnvrt(double x, int ipwr2) {
	return ((int)(ipwr2 * x));
}

/*--------------------------------------------------------------------
c       set ith element of sparse vector (v, iv) with
c       nzv nonzeros to val
c-------------------------------------------------------------------*/
static void vecset(
	int n,
	double v[],	/* v[1:*] */
	int iv[],	/* iv[1:*] */
	int *nzv,
	int i,
	double val)
{
	int k;
	boolean set;

	set = FALSE;
	for (k = 1; k <= *nzv; k++) {
		if (iv[k] == i) {
			v[k] = val;
			set  = TRUE;
		}
	}
	if (set == FALSE) {
		*nzv = *nzv + 1;
		v[*nzv] = val;
		iv[*nzv] = i;
	}
}
