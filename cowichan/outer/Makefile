include ../make.defs

.PHONY: all
all: outer

outer: outer.o
	$(CXX) -o $@ $? $(LIB) 

outer.o: outer.cpp
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
	rm -f outer
	rm -f outer.o
