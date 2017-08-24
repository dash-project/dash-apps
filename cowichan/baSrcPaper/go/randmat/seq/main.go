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

func randmat(nrows, ncols, s int) {
	LCG_A := uint32(1664525)
	LCG_C := uint32(1013904223)
	for i := 0; i < nrows; i++ {
		seed := uint32(s) + uint32(i)
		for j := 0; j < ncols; j++ {
			seed = LCG_A*seed + LCG_C
			matrix[i][j] = int(seed % 100)
		}
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
