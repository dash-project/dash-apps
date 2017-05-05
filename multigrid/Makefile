include make.defs

.phony: all
all: multigrid heat_equation2d heat_equation3d

multigrid:  multigrid.cpp
	$(CXX) -c $(INC) `libpng-config --cflags` $?
	$(CXX) -o $@ $@.o $(LIB) `libpng-config --ldflags` -lhwloc -lnuma

heat_equation2d:  heat_equation2d.cpp
	$(CXX) -c $(INC) $?
	$(CXX) -o $@ $@.o $(LIB) -lhwloc -lnuma

heat_equation3d:  heat_equation3d.cpp
	$(CXX) -c $(INC) $?
	$(CXX) -o $@ $@.o $(LIB) -lhwloc -lnuma

.phony: printenv
printenv :
	@echo "CXX           = $(CXX)"
	@echo "DART_IMPL     = $(DART_IMPL)"
	@echo "DASH_ROOT     = $(DASH_ROOT)"
	@echo "INC           = $(INC)"
	@echo "LIB           = $(LIB)"

.phony: clean
clean:
	rm -f heat_equation*d multigrid halo_heat_eqn *.o *.png
