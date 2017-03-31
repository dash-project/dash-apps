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
	"runtime"
)

type ByteMatrix struct {
	Rows, Cols int
	array      []byte
}

func NewByteMatrix(r, c int) *ByteMatrix {
	return &ByteMatrix{r, c, make([]byte, r*c)}
}

func (m *ByteMatrix) Row(i int) []byte {
	return m.array[i*m.Cols : (i+1)*m.Cols]
}

const (
	LCG_A = 1664525
	LCG_C = 1013904223
)

var (
	is_bench = flag.Bool("is_bench", false, "")
)

func randmat(nrows, ncols int, s uint32) *ByteMatrix {
	matrix := NewByteMatrix(nrows, ncols)

	work := make(chan int)

	go func() {
		for i := 0; i < nrows; i++ {
			work <- i
		}
		close(work)
	}()

	done := make(chan bool)
	NP := runtime.GOMAXPROCS(0)

	for i := 0; i < NP; i++ {
		go func() {
			for i := range work {
				seed := s + uint32(i)
				row := matrix.Row(i)
				for j := range row {
					seed = LCG_A*seed + LCG_C
					row[j] = byte(seed%100) % 100
				}
			}
			done <- true
		}()
	}

	for i := 0; i < NP; i++ {
		<-done
	}

	return matrix
}

func main() {
	flag.Parse()

	var nrows, ncols int
	var seed uint32

	fmt.Scan(&nrows)
	fmt.Scan(&ncols)
	fmt.Scan(&seed)
	matrix := randmat(nrows, ncols, seed)

	if !*is_bench {
		for i := 0; i < nrows; i++ {
			row := matrix.Row(i)
			for j := range row {
				fmt.Printf("%d ", row[j])
			}
			fmt.Printf("\n")
		}
		fmt.Printf("\n")
	}
}
