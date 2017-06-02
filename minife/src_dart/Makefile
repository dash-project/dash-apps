#-----------------------------------------------------------------------

MINIFE_TYPES =  \
        -DMINIFE_SCALAR=double   \
        -DMINIFE_LOCAL_ORDINAL=int      \
        -DMINIFE_GLOBAL_ORDINAL=int

MINIFE_MATRIX_TYPE = -DMINIFE_CSR_MATRIX
# MINIFE_MATRIX_TYPE = -DMINIFE_ELL_MATRIX

#-----------------------------------------------------------------------

DARTROOT=$(HOME)/opt/dash-0.3.0/

CFLAGS = -O3 -g -DDART_FULL_ALLOC
CXXFLAGS = $(CFLAGS) -std=c++11

CPPFLAGS = -I. -I../utils -I../fem -I$(DARTROOT)/include $(MINIFE_TYPES) $(MINIFE_MATRIX_TYPE) \
	-DHAVE_MPI -DHAVE_DART -DMPICH_IGNORE_CXX_SEEK \
	-DMINIFE_REPORT_RUSAGE

LDFLAGS = -L$(DARTROOT)/lib
LIBS=-ldash-mpi -ldart-mpi -ldart-base -lnuma -lhwloc

CXX=mpicxx 
CC=mpicc

#CXX=g++
#CC=gcc

include make_targets

