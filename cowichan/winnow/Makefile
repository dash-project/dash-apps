#include ../make.defs

.PHONY: all
all: main

winnow: winnow.cpp
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
	-rm main
#-rm winnow
#-rm winnow.o

