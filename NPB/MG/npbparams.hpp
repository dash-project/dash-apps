/* CLASS = B */
/*
c  This file is generated automatically by the setparams utility.
c  It sets the number of processors and the class_npb of the NPB
c  in this directory. Do not modify it by hand.
*/
#define	NX_DEFAULT	256
#define	NY_DEFAULT	256
#define	NZ_DEFAULT	256
#define	NIT_DEFAULT	20
#define	LM	8
#define	LT_DEFAULT	8
#define	DEBUG_DEFAULT	0
#define	NDIM1	8
#define	NDIM2	8
#define	NDIM3	8
#define	CONVERTDOUBLE	FALSE
#define COMPILETIME "10 Oct 2019"
#define NPBVERSION "4.0"
#define CS1 "dash-mpicxx"
#define CS2 "$(CC)"
#define CS3 "-L$(HOME)/opt/dash-0.4.0/lib"
#define CS4 "-I../common -I$(HOME)/opt/dash-0.4.0/include"
#define CS5 "-std=c++14 -O3"
#define CS6 "-lm -ldash-mpi -ldart-mpi -ldart-base -lhwloc -lnuma -lpthread"
#define CS7 "randdp"
