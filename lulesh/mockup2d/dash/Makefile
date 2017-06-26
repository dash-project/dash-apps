
BUILD = manual

CXX = mpicxx

ifeq ($(BUILD),cmake)
  DASH_ROOT = $(HOME)/opt/dash-0.3.0/
  DASH_LIB = -L$(DASH_ROOT)/lib -ldash-mpi -ldart-mpi -ldart-base \
#                                -lhwloc \
#                                -lnuma  \
#                                -lpapi
  DASH_INC = -I$(DASH_ROOT)/include
endif

ifeq ($(BUILD),manual)
  DASH_LIB = $(HOME)/code/dash/dash/src/libdash.a           \
	           $(HOME)/code/dash/dart-impl/mpi/src/libdart.a

  DASH_INC = -I$(HOME)/code/dash/dash/include/              \
             -I$(HOME)/code/dash/dart-if/include/
endif

all : lulesh-mock.o domain.o comm.o
	$(CXX) -o lulesh-mock $^ $(DASH_LIB)

%.o : %.cc
	$(CXX) -c -o $@ $(DASH_INC) $+

clean :
	rm -f ./*.o
	rm -f *~
	rm -f ./lulesh-mock
