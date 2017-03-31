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

type ByteMatrix struct {
	Rows, Cols uint32
	array      []byte
}

func NewByteMatrix(r, c uint32) *ByteMatrix {
	return &ByteMatrix{r, c, make([]byte, r*c)}
}

func WrapBytes(r, c uint32, bytes []byte) *ByteMatrix {
	return &ByteMatrix{r, c, bytes}
}

func (m *ByteMatrix) Row(i uint32) []byte {
	return m.array[i*m.Cols : (i+1)*m.Cols]
}

func (m *ByteMatrix) Bytes() []byte {
	return m.array[0 : m.Rows*m.Cols]
}

var mask [][]bool

func thresh(m *ByteMatrix, nrows, ncols, percent uint32) {
	var hist [100]int

	for _, v := range m.Bytes() {
		hist[v]++
	}

	count := (nrows * ncols * percent) / 100
	prefixsum := 0
	threshold := 99

	for ; threshold > 0; threshold-- {
		prefixsum += hist[threshold]
		if prefixsum > int(count) {
			break
		}
	}
	for i := uint32(0); i < nrows; i++ {
		row := m.Row(i)
		for j := range row {
			mask[i][j] = row[j] >= byte(threshold)
		}
	}
}

func main() {
	var nrows, ncols, percent uint32

	flag.Parse()

	fmt.Scanf("%d%d", &nrows, &ncols)
  mask = make ([][]bool, nrows)
  for i := range mask {
    mask [i] = make ([]bool, ncols)
  }

	m := WrapBytes(nrows, ncols, make([]byte, nrows*ncols))

	if !*is_bench {
		for i := uint32(0); i < nrows; i++ {
			row := m.Row(i)
			for j := range row {
				fmt.Scanf("%d", &row[j])
			}
		}
	}

	fmt.Scanf("%d", &percent)

	thresh(m, nrows, ncols, percent)

	if !*is_bench {
		for i := uint32(0); i < nrows; i++ {
			for j := uint32(0); j < ncols; j++ {
				if mask[i][j] {
					fmt.Printf("1 ")
				} else {
					fmt.Printf("0 ")
				}
			}
			fmt.Printf("\n")
		}
		fmt.Printf("\n")
	}
}
