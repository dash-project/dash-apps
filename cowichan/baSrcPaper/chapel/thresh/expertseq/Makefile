all: main
	./main < main.in > main.out
	diff main.out main.g
	./main < main2.in > main2.out
	diff main2.out main2.g
	./main < main3.in > main3.out
	diff main3.out main3.g
	./main < main4.in > main4.out
	diff main4.out main4.g
	./main < main5.in > main5.out
	diff main5.out main5.g
	./main < main6.in > main6.out
	diff main6.out main6.g
	./main < main.in
	./main --is_bench < main.in

main: main.chpl
	chpl --fast main.chpl -o main
