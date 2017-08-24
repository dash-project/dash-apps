all: main
	./main < main.in
	echo "1 4 8" | ./main --nproc 4
	echo "1 5 8" | ./main --nproc 4
	echo "2 5 8" | ./main --nproc 4
	echo "2 5 9" | ./main --nproc 4
	echo "2 5 9" | ./main --nproc 4 --is_bench

# flags: -cilk-profile, -cilk-span -nproc 4
main: main.cpp
	icpc ${CFLAGS} -O2 -o main $<
