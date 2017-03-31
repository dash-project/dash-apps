/*
 * thresh: histogram thresholding
 *
 * input:
 *   matrix: the integer matrix to be thresholded
 *   nrows, ncols: the number of rows and columns
 *   percent: the percentage of cells to retain
 *
 * output:
 *   mask: a boolean matrix filled with true for cells that are kept
 *
 */
package main

import (
	"flag"
	"fmt"
)

var is_bench = flag.Bool("is_bench", false, "")
var matrix [][]byte
var mask [][]byte
var histogram [100]int

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func thresh(nrows, ncols int, percent int) {
	var nmax int = 0
	for i := 0; i < nrows; i++ {
		for j := 0; j < ncols; j++ {
			nmax = max(nmax, int(matrix[i][j]))
		}
	}

	for i := 0; i < nrows; i++ {
		for j := 0; j < ncols; j++ {
			histogram[matrix[i][j]]++
		}
	}

	var count int = (nrows * ncols * percent) / 100
	var prefixsum int = 0
	var threshold int = nmax

	for i := nmax; i >= 0 && prefixsum <= count; i-- {
		prefixsum += histogram[i]
		threshold = i
	}

	for i := 0; i < nrows; i++ {
		for j := 0; j < ncols; j++ {
			if int(matrix[i][j]) >= threshold {
				mask[i][j] = 1
			}
		}
	}
}

func main() {
	var nrows, ncols, percent int

	flag.Parse()

	fmt.Scanf("%d%d", &nrows, &ncols)
  matrix = make ([][]byte, nrows)
  for i := range matrix {
    matrix [i] = make ([]byte, ncols)
  }

  mask = make ([][]byte, nrows)
  for i := range mask {
    mask [i] = make ([]byte, ncols)
  }


	if !*is_bench {
		for i := 0; i < nrows; i++ {
			for j := 0; j < ncols; j++ {
				fmt.Scanf("%d", &matrix[i][j])
			}
		}
	}

	fmt.Scanf("%d", &percent)

	thresh(nrows, ncols, percent)

	if !*is_bench {
		for i := 0; i < nrows; i++ {
			for j := 0; j < ncols; j++ {
				fmt.Printf("%d ", mask[i][j])
			}
			fmt.Printf("\n")
		}
		fmt.Printf("\n")
	}
}
