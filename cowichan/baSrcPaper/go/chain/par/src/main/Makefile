all: main
	./main < main.in
	./main < main.in --is_bench

DEPS=../all/randmat.go ../all/thresh.go ../all/winnow.go ../all/outer.go \
		../all/product.go
main: main.go $(DEPS)
	GOPATH="`pwd`/../.." go build
