all: main
	./$< < $<.in > main.out --threads 4
	diff main.out main.gold
	./$< < $<.in --threads 4
	./$< < $<.in --threads 4 --is_bench

main: randmat.o thresh.o winnow.o outer.o product.o

CXXFLAGS+=-std=c++0x -O2 -m64 -ltbb -lpthread -lrt
LDFLAGS+=-ltbb
CXX=icpc
