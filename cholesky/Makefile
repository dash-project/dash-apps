

CXX = dash-mpicxx

CXXFLAGS = -std=c++14

LDFLAGS = -llapacke -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl -ldart-tasking

OBJ = cholesky cholesky_tasks cholesky_tasks_prefetch

HEADER = common.h MatrixBlock.h

all: $(OBJ)

cholesky: main.cpp $(HEADER)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

cholesky_tasks: main_tasks.cpp $(HEADER)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

cholesky_tasks_prefetch: main_tasks_prefetch.cpp $(HEADER)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm *.o $(OBJ)
