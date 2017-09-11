all: main
	./main < main.in > main.out
	diff main.out main.gold
	./main < main.in
	./main < main.in --is_bench

main: main.6
	go tool 6l -o main main.6

main.6: main.go
	go tool 6g -O2 main.go
