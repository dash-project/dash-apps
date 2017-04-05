include ../make.defs

.PHONY: all
all: thresh

thresh: thresh.cpp
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
	-rm thresh
	-rm thresh.o
