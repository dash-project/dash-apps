CXXFLAGS+=-std=c++0x -O2 -m64 -ltbb -lpthread -lrt -lrt
all: main
	echo "10 9 8" | ./main --threads 4
	echo "1 4 8" | ./main --threads 4
	echo "1 5 8" | ./main --threads 4
	echo "2 5 8" | ./main --threads 4
	echo "2 5 9" | ./main --threads 4
	echo "2 5 9" | ./main --is_bench --threads 4

CXXFLAGS+=-std=c++0x -O2
LDFLAGS+=-ltbb
CXX=icpc
