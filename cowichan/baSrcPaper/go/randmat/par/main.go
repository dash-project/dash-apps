/*
 * randmat: random number generation
 *
 * input:
 *   nrows, ncols: the number of rows and columns
 *   s: the seed
 *
 * output:
 *   martix: a nrows x ncols integer matrix
 *
 */
package main

import (
	"flag"
	"fmt"
)

var is_bench = flag.Bool("is_bench", false, "")
var matrix [][]int

func randvec(row, n, seed int, done chan bool) {
	LCG_A := uint32(1664525)
	LCG_C := uint32(1013904223)
  s := uint32(seed)
  for j := 0; j < n; j++ {
	  s = LCG_A*s + LCG_C
		matrix[row][j] = int(s % 100)
	}
	done <- true
}

func parallel_for(begin, end, ncols, s int, done chan bool) {
	if begin+1 == end {
		randvec(begin, ncols, s+begin, done)
	} else {
		middle := begin + (end-begin)/2
		go parallel_for(begin, middle, ncols, s, done)
		parallel_for(middle, end, ncols, s, done)
	}
}

func randmat(nrows, ncols, s int) {
	done := make(chan bool)
	// parallel for on rows
	go parallel_for(0, nrows, ncols, s, done)
	for i := 0; i < nrows; i++ {
		<-done
	}
}

func main() {
	var nrows, ncols, s int

	flag.Parse()

	fmt.Scanf("%d%d%d", &nrows, &ncols, &s)
	matrix = make([][]int, nrows)
	for i := range matrix {
		matrix[i] = make([]int, ncols)
	}

	randmat(nrows, ncols, s)

	if !*is_bench {
		for i := 0; i < nrows; i++ {
			for j := 0; j < ncols; j++ {
				fmt.Printf("%d ", matrix[i][j])
			}
			fmt.Printf("\n")
		}
		fmt.Printf("\n")
	}
}
