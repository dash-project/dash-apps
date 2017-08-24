CXXFLAGS+=-std=c++0x -O2 -m64 -ltbb -lpthread -lrt -lrt
all: main
	echo "10 9 8" | ./main
	echo "1 4 8" | ./main
	echo "1 5 8" | ./main
	echo "2 5 8" | ./main
	echo "2 5 9" | ./main
	echo "2 5 9" | ./main --is_bench

CXXFLAGS+=-std=c++0x -O2
LDFLAGS+=-ltbb
CXX=icpc
