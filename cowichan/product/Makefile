include ../make.defs

.PHONY: all
all: product

product: product.cpp product.h
	$(CXX) -c $(INC) $<
	$(CXX) -o $@ $@.o $(LIB) 


.PHONY: printenv
printenv :
	@echo "CXX           = $(CXX)"
	@echo "DART_IMPL     = $(DART_IMPL)"
	@echo "DASH_ROOT     = $(DASH_ROOT)"
	@echo "INC           = $(INC)"
	@echo "LIB           = $(LIB)"

.PHONY: clean
clean:
	-rm product
	-rm product.o
