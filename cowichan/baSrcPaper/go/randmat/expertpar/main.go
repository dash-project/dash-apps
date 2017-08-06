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
  "bufio"
  "os"
)

// #include <time.h>
// #include <stdio.h>
import "C"

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
  
  var start, stop C.struct_timespec
  var accum C.double
  
  if( C.clock_gettime( C.CLOCK_MONOTONIC_RAW, &start) == -1 ) {
    C.perror( C.CString("clock gettime error 1") );
    return
  }
  
	matrix := randmat(nrows, ncols, seed)

  if( C.clock_gettime( C.CLOCK_MONOTONIC_RAW, &stop) == -1 ) {
    C.perror( C.CString("clock gettime error 1") );
    return
  }
  
  accum = C.double( C.long(stop.tv_sec) - C.long(start.tv_sec) ) + C.double(( C.long(stop.tv_nsec) - C.long(start.tv_nsec))) / C.double(1e9);
  
  file, err := os.OpenFile("./measurements.txt", os.O_WRONLY|os.O_CREATE|os.O_APPEND, 0666)
  
  if err != nil {
      fmt.Println("File does not exists or cannot be created")
      os.Exit(1)
  }
  
  w := bufio.NewWriter(file)
  
  // Lang, Problem, rows, cols, thresh, winnow_nelts, jobs, time
  fmt.Fprintf(w, "Go,Randmat,%d, %d, , , %d,%.9f,isBench:%t\n", nrows, ncols, runtime.GOMAXPROCS(0), accum, *is_bench )
  
  w.Flush()
  file.Close()
  
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