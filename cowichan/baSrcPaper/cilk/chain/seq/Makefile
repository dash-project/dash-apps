all: main
	./main < main.in > main.out
	diff main.out main.gold
	./main --nproc 4 < main.in
	./main --is_bench < main.in

# flags: --cpp-profile, --cpp-span
main: randmat.cpp thresh.cpp winnow.cpp outer.cpp product.cpp main.cpp
	icpc ${CFLAGS} -O2 -o main $^
