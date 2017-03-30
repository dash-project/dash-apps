include ../make.defs

.PHONY: all
all: winnow

winnow: winnow.o
	$(CXX) -o $@ $? $(LIB) 

winnow.o: winnow.cpp
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
	rm -f winnow
	rm -f winnow.o
