/*
 * product: matrix-vector product
 *
 * input:
 *   nelts: the number of elements
 *   matrix: the real matrix
 *   vector: the real vector
 *
 * output:
 *   result: a real vector, whose values are the result of the product
 */
package main

import (
	"flag"
	"fmt"
)

var is_bench = flag.Bool("is_bench", false, "")

type double float64

var matrix [][]double
var vector []double
var result []double

func fill_result_impl(begin, end, ncols int, done chan bool) {
	if begin+1 == end {
		var sum double = 0
		for j := 0; j < ncols; j++ {
			sum += matrix[begin][j] * vector[j]
		}
		result[begin] = sum
		done <- true
	} else {
		middle := begin + (end-begin)/2
		go fill_result_impl(begin, middle, ncols, done)
		fill_result_impl(middle, end, ncols, done)
	}
}

func fill_result(nrows, ncols int) {
	done := make(chan bool)
	// parallel for on rows
	go fill_result_impl(0, nrows, ncols, done)
	for i := 0; i < nrows; i++ {
		<-done
	}
}

func product(nelts int) {
	fill_result(nelts, nelts)
}

func read_integer() int {
	var value int
	for true {
		var read, _ = fmt.Scanf("%d", &value)
		if read == 1 {
			break
		}
	}
	return value
}

func read_double() double {
	var value double
	for true {
		var read, _ = fmt.Scanf("%g", &value)
		if read == 1 {
			break
		}
	}
	return value
}

func read_matrix(nelts int) {
	for i := 0; i < nelts; i++ {
		for j := 0; j < nelts; j++ {
			matrix[i][j] = read_double()
		}
	}
}

func read_vector(nelts int) {
	for i := 0; i < nelts; i++ {
		vector[i] = read_double()
	}
}

func main() {
	var nelts int

	flag.Parse()

	nelts = read_integer()
  matrix = make ([][]double, nelts)
  for i := range matrix {
    matrix [i] = make ([]double, nelts)
  }
  vector = make ([]double, nelts)
  result = make ([]double, nelts)


	if !*is_bench {
		read_matrix(nelts)
		read_vector(nelts)
	}

	product(nelts)

	if !*is_bench {
		fmt.Printf("%d\n", nelts)
		for i := 0; i < nelts; i++ {
			fmt.Printf("%g ", result[i])
		}
		fmt.Printf("\n")
	}
}
