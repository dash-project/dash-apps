

CXX = dash-mpicxx

#CXXFLAGS = -std=c++14 -mkl=sequential -O3
CXXFLAGS = -std=c++14 -ggdb

#LDFLAGS = -mkl=sequential -lpthread -lm -ldl
LDFLAGS = -llapacke -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl
DART_TASK_LIB = -ldart-tasking

EXTRAE_ROOT= $(HOME)/opt/extrae-3.5.2/
EXTRAE_INC = $(EXTRA_ROOT)/include/
EXTRAE_LIB = $(EXTRA_ROOT)/lib

OBJ = cholesky cholesky_prefetch cholesky_prefetch_omp cholesky_prefetch_extrae cholesky_tasks cholesky_tasks_prefetch cholesky_tasks_prefetch_extrae cholesky_tasks_copyin cholesky_tasks_copyin_extrae
MAIN = main.cpp

HEADER = common.h MatrixBlock.h

.PHONY: clean

all: $(OBJ)

cholesky: $(MAIN) $(HEADER) Cholesky.h Makefile
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

cholesky_prefetch: $(MAIN) $(HEADER) CholeskyPrefetch.h Makefile
	$(CXX) $(CXXFLAGS) -DDASH_PREFETCH -o $@ $< $(LDFLAGS)

cholesky_prefetch_extrae: $(MAIN) $(HEADER) CholeskyPrefetch.h
	$(CXX) $(CXXFLAGS) -DDASH_PREFETCH -DUSE_EXTRAE -I$(EXTRAE_INC) -o $@ $< $(LDFLAGS) $(DART_TASK_LIB) -L$(EXTRAE_LIB) -lptmpitrace

cholesky_prefetch_omp: $(MAIN) $(HEADER) CholeskyOMPPrefetch.h Makefile
	$(CXX) $(CXXFLAGS) -DDASH_OMP_PREFETCH -fopenmp -o $@ $< $(LDFLAGS)

cholesky_tasks: $(MAIN) $(HEADER) CholeskyTasks.h Makefile
	$(CXX) $(CXXFLAGS) -DDASH_TASKS -o $@ $< $(LDFLAGS) $(DART_TASK_LIB)

cholesky_tasks_prefetch: $(MAIN) $(HEADER) CholeskyTasksPrefetch.h Makefile
	$(CXX) $(CXXFLAGS) -DDASH_TASKS_PREFETCH -o $@ $< $(LDFLAGS) $(DART_TASK_LIB)

cholesky_tasks_prefetch_extrae: $(MAIN) $(HEADER) CholeskyTasksPrefetch.h
	$(CXX) $(CXXFLAGS) -DDASH_TASKS_PREFETCH -DUSE_EXTRAE -I$(EXTRAE_INC) -o $@ $< $(LDFLAGS) $(DART_TASK_LIB) -L$(EXTRAE_LIB) -lptmpitrace

cholesky_tasks_copyin: $(MAIN) $(HEADER) CholeskyTasksCopyIn.h Makefile
	$(CXX) $(CXXFLAGS) -DDASH_TASKS_COPYIN -o $@ $< $(LDFLAGS) $(DART_TASK_LIB)

cholesky_tasks_copyin_extrae: $(MAIN) $(HEADER) CholeskyTasksCopyIn.h
	$(CXX) $(CXXFLAGS) -DDASH_TASKS_COPYIN -DUSE_EXTRAE -I$(EXTRAE_INC) -o $@ $< $(LDFLAGS) $(DART_TASK_LIB) -L$(EXTRAE_LIB) -lptmpitrace

clean:
	rm -rf $(OBJ)
