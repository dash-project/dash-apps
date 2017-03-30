include ../make.defs

.PHONY: all
all: randmat
#mpirun -n 4 ./$@ < main.in > rand.out2
#diff -ws rand.out2 main.gold

randmat:  randmat.o
	$(CXX) -o $@ $? $(LIB)

randmat.o:randmat.cpp 
	$(CXX) -c $(INC) $?


.PHONY: printenv
printenv :
	@echo "CXX           = $(CXX)"
	@echo "DART_IMPL     = $(DART_IMPL)"
	@echo "DASH_ROOT     = $(DASH_ROOT)"
	@echo "INC           = $(INC)"
	@echo "LIB           = $(LIB)"

.PHONY: clean
clean:
	rm -f randmat
	rm -f randmat.o
