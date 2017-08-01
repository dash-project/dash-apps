
BUILD = manual

CXX = mpicxx

#ifeq ($(BUILD),cmake)
  LIB = -L$(DASH_LIB) -L$(DART_LIB) -ldash-mpi -ldart-mpi -ldart-base \
#                                -lhwloc \
#                                -lnuma  \
#                                -lpapi
  INC = -I$(DASH_INC)
#endif

#ifeq ($(BUILD),manual)
#  DASH_LIB = $(HOME)/code/dash/dash/src/libdash.a           \
	           $(HOME)/code/dash/dart-impl/mpi/src/libdart.a

#  DASH_INC = -I$(HOME)/code/dash/dash/include/              \
             -I$(HOME)/code/dash/dart-if/include/
#endif

all : lulesh-mock.o domain.o
	$(CXX) -o lulesh-mock $^ $(LIB)

%.o : %.cc
	$(CXX) -c -std=c++11 -o $@ $(INC) $+

clean :
	rm -f ./*.o
	rm -f *~
	rm -f ./lulesh-mock
