CXXFLAGS+=-std=c++0x -O2 -m64 -ltbb -lpthread -lrt -lrt
all: main
	./$< < $<.in > main.out
	diff main.out main.gold
	./$< < $<2.in > main2.out
	diff main2.out main2.gold
	./$< < $<3.in > main3.out
	diff main3.out main3.gold
	./$< < $<4.in > main4.out
	diff main4.out main4.gold
	./$< < $<.in > main.out
	./$< < $<.in
	./$< < $<.in --is_bench

CXXFLAGS+=-std=c++0x -O2
LDFLAGS+=-ltbb
CXX=icpc
