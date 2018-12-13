

CXX = dash-mpicxx

CXXFLAGS = -std=c++14 -mkl=sequential -O3 -pthread -mt_mpi
CXXFLAGS = -std=c++14 -ggdb -I${MKLROOT}/include/

LDFLAGS = -mkl=sequential -pthread -lm -ldl -mt_mpi
LDFLAGS = -llapacke -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl
DART_TASK_LIB = -ldart-tasking

COREBLAS_DIR=$(HOME)/opt/plasma-2.8.0/

CXXFLAGS += -I$(COREBLAS_DIR)/include
LDFLAGS  += -L$(COREBLAS_DIR)/lib -lcoreblas

EXTRAE_ROOT = /work/xg17i034/x10288/opt/extrae-3.5.4/
EXTRAE_ROOT= $(HOME)/opt/extrae-3.5.4/
EXTRAE_INC = $(EXTRAE_ROOT)/include/
EXTRAE_LIB = $(EXTRAE_ROOT)/lib

OBJ = qrfact_tasks_copyin

MAIN = main.cpp

HEADER = common.h MatrixBlock.h

.PHONY: clean

all: $(OBJ)

qrfact_tasks_copyin: $(MAIN) $(HEADER) QrFactTasksCopyin.h Makefile
	$(CXX) $(CXXFLAGS) -DDASH_TASKS_COPYIN -o $@ $< $(LDFLAGS) $(DART_TASK_LIB)

clean:
	rm -rf $(OBJ)
