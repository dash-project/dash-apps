all: main
	./main < main.in > main.out
	diff main.out main.gold
	./main < main2.in > main2.out
	diff main2.out main2.gold
	./main < main3.in > main3.out
	diff main3.out main3.gold
	./main < main4.in > main4.out
	diff main4.out main4.gold
	./main < main.in
	./main < main.in --is_bench

main: main.go
	go build -o main
