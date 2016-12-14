include ../make.defs

.phony: all
all: astro heic1620a.tif

.phon<: diff
diff:  astro-assignment.cpp astro-solution.cpp
	diff  astro-assignment.cpp astro-solution.cpp

heic1620a.tif:
	@echo "  --> please download laegest version of heic1620a.tif at http://www.spacetelescope.org/images/heic1620a/"

astro:  astro.cpp misc.h
	$(CXX) -c $(INC) -DWITHSDL $?
	$(CXX) -o $@ $@.o $(LIB) -lSDL -ltiff


.phony: printenv
printenv :
	@echo "CXX           = $(CXX)"
	@echo "DART_IMPL     = $(DART_IMPL)"
	@echo "DASH_ROOT     = $(DASH_ROOT)"
	@echo "INC           = $(INC)"
	@echo "LIB           = $(LIB)"

.phony: clean
clean:
	rm -f astro *.o
