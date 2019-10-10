#define CLASS 'B'
/*
   This file is generated automatically by the setparams utility.
   It sets the number of processors and the class of the NPB
   in this directory. Do not modify it by hand.   */
   
#define COMPILETIME "02 Oct 2019"
#define NPBVERSION "4.0"
#define CC "mpicxx"
#define CFLAGS "-std=c++14 -O3"
#define CLINK "$(CC)"
#define CLINKFLAGS "-lm -ldash-mpi -ldart-mpi -ldart-base -lhwloc -lnuma -lpthread"
#define C_LIB "-L$(HOME)/opt/dash-0.4.0/lib"
#define C_INC "-I../common -I$(HOME)/opt/dash-0.4.0/include"
