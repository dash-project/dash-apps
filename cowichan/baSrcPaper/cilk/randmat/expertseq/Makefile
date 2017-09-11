all: main
	./main < main.in
	echo "1 4 8" | ./main
	echo "1 5 8" | ./main
	echo "2 5 8" | ./main
	echo "2 5 9" | ./main
	echo "2 5 9" | ./main --is_bench

main: main.cpp
	icpc ${CFLAGS} -O2 -o main main.cpp
