include ../make.defs

.phony: all
all: hello hello_array hello_minimum

hello:  hello.cpp
	$(CXX) -c $(INC) $?
	$(CXX) -o $@ $@.o $(LIB) -lhwloc -lnuma

hello_array:  hello_array.cpp
	$(CXX) -c $(INC) $?
	$(CXX) -o $@ $@.o $(LIB) -lhwloc -lnuma

hello_minimum:  hello_minimum.cpp
	$(CXX) -c $(INC) $?
	$(CXX) -o $@ $@.o $(LIB) -lhwloc -lnuma

.phony: printenv
printenv :
	@echo "CXX           = $(CXX)"
	@echo "DART_IMPL     = $(DART_IMPL)"
	@echo "DASH_ROOT     = $(DASH_ROOT)"
	@echo "INC           = $(INC)"
	@echo "LIB           = $(LIB)"

.phony: clean
clean:
	rm -f hello hello_array hello_minimum *.o
