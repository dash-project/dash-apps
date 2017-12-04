

CXX = dash-mpicxx

CXXFLAGS = -std=c++14 -mkl=sequential -O3

LDFLAGS = -mkl=sequential -lpthread -lm -ldl
DART_TASK_LIB = -ldart-tasking

EXTRAE_INC = /work/xg17i034/x10288/opt/extrae-3.5.2/include/

OBJ = cholesky cholesky_tasks cholesky_tasks_prefetch cholesky_tasks_prefetch_extrae
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

cholesky_tasks_prefetch_extrae: $(MAIN) $(HEADER) CholeskyTasksPrefetch.h
	$(CXX) $(CXXFLAGS) -DDASH_TASKS_PREFETCH -DUSE_EXTRAE -I$(EXTRAE_INC) -L/work/xg17i034/x10288/opt/extrae-3.5.2/lib -lptmpitrace -o $@ $< $(LDFLAGS) $(DART_TASK_LIB)
clean:
	rm -rf $(OBJ)
