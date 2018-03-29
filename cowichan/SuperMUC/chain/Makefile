include ../make.defs

.PHONY: all
all: chain
#mpirun -n 4 ./$@ < main.in > rand.out2
#diff -ws rand.out2 main.gold

chain: chain.cpp ../randmat/randmat.h ../thresh/thresh.h ../winnow/winnow.h ../outer/outer.h ../product/product.h
	$(CXX) -c $(INC) $<
	$(CXX) -o $@ $@.o $(LIB)


.PHONY: printenv
printenv:
	@echo "CXX           = $(CXX)"
	@echo "DART_IMPL     = $(DART_IMPL)"
	@echo "DASH_ROOT     = $(DASH_ROOT)"
	@echo "INC           = $(INC)"
	@echo "LIB           = $(LIB)"

.PHONY: clean
clean:
	-rm chain
	-rm chain.o
