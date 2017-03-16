include ../make.defs

.phony: all
all: winnow

winnow: winnow.cpp
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
	rm -f winnow
	rm -f winnow.o
