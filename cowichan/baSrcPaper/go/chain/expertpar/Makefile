all: main
	./main < main.in > main.out; 
	diff -Nuar main.out main.gold

main: main.go
	go build -o main
