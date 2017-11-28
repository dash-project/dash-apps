

CXX = dash-mpicxx

CXXFLAGS = -std=c++14

LDFLAGS = -llapacke -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl
DART_TASK_LIB = -ldart-tasking

OBJ = cholesky cholesky_tasks cholesky_tasks_prefetch
MAIN = main.cpp

HEADER = common.h MatrixBlock.h

.PHONY: clean

all: $(OBJ)

cholesky: $(MAIN) $(HEADER) Cholesky.h
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

cholesky_tasks: $(MAIN) $(HEADER) CholeskyTasks.h
	$(CXX) $(CXXFLAGS) -DDASH_TASKS -o $@ $< $(LDFLAGS) $(DART_TASK_LIB)

cholesky_tasks_prefetch: $(MAIN) $(HEADER) CholeskyTasksPrefetch.h
	$(CXX) $(CXXFLAGS) -DDASH_TASKS_PREFETCH -o $@ $< $(LDFLAGS) $(DART_TASK_LIB)

clean:
	rm -rf $(OBJ)
