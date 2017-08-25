all: main
	./main < main.in > main.out
	diff main.out main.gold
	./main --nproc 4 < main.in
	./main --nproc 4 --is_bench < main.in

.cpp.o:
	icpc ${CFLAGS} -O2 -c $<

main: main.o randmat.o thresh.o winnow.o outer.o product.o
	icpc ${CFLAGS} -O2 -o main $^
