CXXFLAGS+=-std=c++0x -O2 -m64 -ltbb -lpthread -lrt -lrt
all: main
	./$< < $<.in --threads 4 > main.out
	diff main.out main.gold
	./$< < $<2.in --threads 4 > main2.out
	diff main2.out main2.gold
	./$< < $<3.in --threads 4 > main3.out
	diff main3.out main3.gold
	./$< < $<4.in --threads 4 > main4.out
	diff main4.out main4.gold
	./$< < $<5.in --threads 4 > main5.out
	diff main5.out main5.gold
	./$< < $<6.in --threads 4 > main6.out
	diff main6.out main6.gold
	./$< < $<.in
	./$< < $<.in --is_bench
CXXFLAGS+=-std=c++0x -O2
LDFLAGS+=-ltbb
CXX=icpc
