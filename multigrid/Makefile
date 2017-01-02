include make.defs

.phony: all
all: multigrid

multigrid:  multigrid.cpp
	$(CXX) -c $(INC) `libpng-config --cflags` $?
	$(CXX) -o $@ $@.o $(LIB) `libpng-config --ldflags`

.phony: printenv
printenv :
	@echo "CXX           = $(CXX)"
	@echo "DART_IMPL     = $(DART_IMPL)"
	@echo "DASH_ROOT     = $(DASH_ROOT)"
	@echo "INC           = $(INC)"
	@echo "LIB           = $(LIB)"

.phony: clean
clean:
	rm -f multigrid halo_heat_eqn *.o *.png
