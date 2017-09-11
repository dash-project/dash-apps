all: main
	./main --dataParTasksPerLocale=4 < main.in
	./main --dataParTasksPerLocale=4 --is_bench < main.in

main: main.chpl
	chpl --fast main.chpl -o main
