all: main
	./main < main.in
	./main --is_bench < main.in

main: main.chpl config.chpl randmat.chpl thresh.chpl winnow.chpl outer.chpl product.chpl
	chpl --fast $^ -o main
