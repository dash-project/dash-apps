all: main
	echo "10 9 8" | ./main
	echo "1 4 6" | ./main
	echo "1 5 6" | ./main
	echo "2 5 6" | ./main
	echo "2 5 9" | ./main
	echo "10 9 6" | ./main --is_bench

main: main.go
	go build -o main
