all: main
	./main < main.in > main.out
	diff main.out main.gold
	./main --nproc 4 < main.in
	./main --is_bench < main.in

main: main.cpp
	icpc ${CFLAGS} -O2 -o main main.cpp
