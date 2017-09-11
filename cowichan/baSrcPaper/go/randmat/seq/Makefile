all: main
	echo "10 9 8" | ./main
	echo "1 4 6" | ./main
	echo "1 5 6" | ./main
	echo "2 5 6" | ./main
	echo "2 5 9" | ./main
	echo "10 9 6" | ./main --is_bench

main: main.6
	go tool 6l -o main main.6

main.6: main.go
	go tool 6g -O2 main.go
