include ../make.defs

.phony: all
all: astro astro-benchmark bug heic1620a.tif

.phon<: diff
diff:  astro-assignment.cpp astro-solution.cpp
	diff  astro-assignment.cpp astro-solution.cpp

heic1620a.tif:
	@echo "  --> please download laegest version of heic1620a.tif at http://www.spacetelescope.org/images/heic1620a/"

astro:  astro.cpp misc.h
	$(CXX) -c $(INC) -DWITHOUTSDL astro.cpp
	$(CXX) -o $@ $@.o $(LIB) -lSDL -ltiff -lhwloc -lnuma

bug:  bug.cpp
	$(CXX) -c $(INC) -DWITHOUTSDL $?
	$(CXX) -o $@ $@.o $(LIB) -lSDL -ltiff -lhwloc -lnuma

astro-benchmark:  astro-benchmark.cpp misc.h
	$(CXX) -c $(INC) -DWITHOUTSDL astro-benchmark.cpp
	$(CXX) -o $@ $@.o $(LIB) -ltiff -lhwloc -lnuma

	
.phony: printenv
printenv :
	@echo "CXX           = $(CXX)"
	@echo "DART_IMPL     = $(DART_IMPL)"
	@echo "DASH_ROOT     = $(DASH_ROOT)"
	@echo "INC           = $(INC)"
	@echo "LIB           = $(LIB)"

.phony: clean
clean:
	rm -f astro astro-benchmark  bug *.gch *.o
