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
	"runtime"
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
	NP := runtime.GOMAXPROCS(0)

	hist_work := make(chan uint32)
	hist_parts := make(chan []int)
	go func() {
		for i := uint32(0); i < nrows; i++ {
			hist_work <- i
		}
		close(hist_work)
	}()

	for i := 0; i < NP; i++ {
		go func() {
			my_hist := make([]int, 100)
			for i := range hist_work {
				row := m.Row(i)
				for j := range row {
					my_hist[row[j]]++
				}
			}
			hist_parts <- my_hist
		}()
	}

	var hist [100]int

	for i := 0; i < NP; i++ {
		my_hist := <-hist_parts
		for j := range my_hist {
			hist[j] += my_hist[j]
		}
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

	mask_work := make(chan uint32)

	go func() {
		for i := uint32(0); i < nrows; i++ {
			mask_work <- i
		}
		close(mask_work)
	}()

	mask_done := make(chan bool)
	for i := 0; i < NP; i++ {
		go func() {
			for i := range mask_work {
				row := m.Row(i)
				for j := range row {
					mask[i][j] = row[j] >= byte(threshold)
				}
			}
			mask_done <- true
		}()
	}

	for i := 0; i < NP; i++ {
		<-mask_done
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


	m := WrapBytes(nrows, ncols, make([]byte, ncols*nrows))

	if !*is_bench {
		for i := uint32(0); i < nrows; i++ {
			row := m.Row(i)
			for j := range row {
				fmt.Scanf("%d", &row[j])
			}
		   fmt.Scanf(" \n")
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
