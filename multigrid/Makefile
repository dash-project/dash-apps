include make.defs

PROG=multigrid3d

.phony: all
all: ${PROG}

.phony: ${PROG}_csv

${PROG}_csv: multigrid3d.cpp allreduce.h minimonitoring.h
	$(CXX) -march=native -DWITHCSVOUTPUT -o $@.o -c $(INC) $<
	$(CXX) -march=native -o $@ $@.o $(LIB) -lrt -lnuma

${PROG}: multigrid3d.cpp allreduce.h minimonitoring.h
	$(CXX) -march=native -c $(INC) $<
	$(CXX) -march=native -o $@ $@.o $(LIB) -lrt -lnuma -lhwloc

multigrid3d_plain.cpp: multigrid3d.cpp
	grep -v -i "minimon"   multigrid3d.cpp > multigrid3d_plain.cpp

multigrid3d_plain: multigrid3d_plain.cpp
	$(CXX) -march=native -c $(INC) $?
	$(CXX) -march=native -o $@ $@.o $(LIB) -lrt -lnuma

multigrid3d_elastic: multigrid3d_elastic.cpp minimonitoring.h
	$(CXX) -march=native -c $(INC) $?
	$(CXX) -march=native -o $@ $@.o $(LIB) -lrt -lnuma

.phony: printenv
printenv :
	@echo "CXX           = $(CXX)"
	@echo "DART_IMPL     = $(DART_IMPL)"
	@echo "DASH_ROOT     = $(DASH_ROOT)"
	@echo "INC           = $(INC)"
	@echo "LIB           = $(LIB)"

.phony: clean
clean:
	rm -f heat_equation*d multigrid multigrid*d multigrid*d+minimon multigrid3d_elastic halo_heat_eqn *.o *.gch
