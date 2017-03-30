include ../make.defs

.phony: all
all: product

product: product.o
	$(CXX) -o $@ $? $(LIB) 

product.o: product.cpp
	$(CXX) -c $(INC) $?


.phony: printenv
printenv :
	@echo "CXX           = $(CXX)"
	@echo "DART_IMPL     = $(DART_IMPL)"
	@echo "DASH_ROOT     = $(DASH_ROOT)"
	@echo "INC           = $(INC)"
	@echo "LIB           = $(LIB)"

.PHONY: clean
clean:
	rm -f product
	rm -f product.o
