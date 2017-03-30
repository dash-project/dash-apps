include ../make.defs

.phony: all
all: thresh

thresh: thresh.o
	$(CXX) -o $@ $? $(LIB) 

thresh.o: thresh.cpp
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
	rm -f thresh
	rm -f thresh.o
