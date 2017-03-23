include ../make.defs

.phony: all
all: randmat

randmat:  randmat.cpp
	$(CXX) -c $(INC) $?
	$(CXX) -o $@ $@.o $(LIB) 


.phony: printenv
printenv :
	@echo "CXX           = $(CXX)"
	@echo "DART_IMPL     = $(DART_IMPL)"
	@echo "DASH_ROOT     = $(DASH_ROOT)"
	@echo "INC           = $(INC)"
	@echo "LIB           = $(LIB)"

.phony: clean
clean:
	rm -f randmat
	rm -f randmat.o
