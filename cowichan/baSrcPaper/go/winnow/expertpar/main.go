/*
 * winnow: weighted point selection
 *
 * input:
 *   matrix: an integer matrix, whose values are used as masses
 *   mask: a boolean matrix showing which points are eligible for
 *     consideration
 *   nrows, ncols: the number of rows and columns
 *   nelts: the number of points to select
 *
 * output:
 *   points: a vector of (x, y) points
 *
 */
package main

import (
	"flag"
	"fmt"
	"runtime"
	"sort"
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

func WrapBytes(r, c int, bytes []byte) *ByteMatrix {
	return &ByteMatrix{r, c, bytes}
}

func NewByteMatrix(r, c int) *ByteMatrix {
	return &ByteMatrix{r, c, make([]byte, r*c)}
}

func (m *ByteMatrix) Row(i int) []byte {
	return m.array[i*m.Cols : (i+1)*m.Cols]
}

func (m *ByteMatrix) Bytes() []byte {
	return m.array[0 : m.Rows*m.Cols]
}

var is_bench = flag.Bool("is_bench", false, "")
var matrix []byte
var mask [][]bool
var points []int

type WinnowPoints struct {
	m *ByteMatrix
	e []int // indexes into the ByteMatrix 'm'
}

func (p *WinnowPoints) Len() int {
	return len(p.e)
}

func (p *WinnowPoints) Swap(i, j int) {
	p.e[i], p.e[j] = p.e[j], p.e[i]
}

func (p *WinnowPoints) Less(i, j int) bool {
	return ArrayLess(p.m.array, p.e[i], p.e[j])
}

func ArrayLess(array []byte, x, y int) bool {
	if array[x] != array[y] {
		return array[x] < array[y]
	}
	return x < y
}

func min(a, b int) int {
    if a < b {
        return a
    }
    return b
}

func WinnowMerge(points chan WinnowPoints) {
	var merged WinnowPoints
	x := <-points
	y := <-points

	new_size := len(x.e) + len(y.e)

	merged.m = x.m
	merged.e = make([]int, new_size)

	j := 0
	k := 0
	for i := 0; i < new_size; i++ {
		if j < len(x.e) && k < len(y.e) {
			if ArrayLess(merged.m.array, x.e[j], y.e[k]) {
				merged.e[i] = x.e[j]
				j++
			} else {
				merged.e[i] = y.e[k]
				k++
			}
		} else if j < len(x.e) {
			merged.e[i] = x.e[j]
			j++
		} else if k < len(y.e) {
			merged.e[i] = y.e[k]
			k++
		}
	}
	points <- merged
}

func Winnow(m *ByteMatrix, nrows, ncols, winnow_nelts int) {
	NP := runtime.GOMAXPROCS(0)
	var values WinnowPoints
	values.m = m

	values_work := make(chan int, 1024)
	values_done := make(chan WinnowPoints, NP)
	values_done <- WinnowPoints{m, make([]int, 0)}

	go func() {
		for i := 0; i < nrows; i++ {
			values_work <- i
		}
		close(values_work)
	}()

	merged := make(chan bool, NP)

	for i := 0; i < NP; i++ {
		go func() {
			WinnowMerge(values_done)
			merged <- true
		}()
	}

	for i := 0; i < NP; i++ {
		go func() {
			var local_indexes []int
			for i := range values_work {
				for j := 0; j < ncols; j++ {
					idx := i*ncols + j
					// if *is_bench {
						// mask[i][j] = ((i * j) % (ncols + 1)) == 1
					// }
					if mask[i][j] {
						local_indexes = append(local_indexes, idx)
					}
				}
			}
			var local_values WinnowPoints
			local_values.m = m
			local_values.e = local_indexes

			sort.Sort(&local_values)
			values_done <- local_values
		}()
	}

	for i := 0; i < NP; i++ {
		<-merged
	}

	values = <-values_done
  
  // fmt.Printf("len:%d\n", values.Len())

	chunk := values.Len() / winnow_nelts

	points = make([]int, winnow_nelts)
	point_work := make(chan int, 1024)
	point_done := make(chan bool)
	go func() {
		for i := 0; i < winnow_nelts; i++ {
			point_work <- i
		}
		close(point_work)
	}()

	for i := 0; i < NP; i++ {
		go func() {
			for i := range point_work {
				points[i] = values.e[i*chunk]
			}
			point_done <- true
		}()
	}

	for i := 0; i < NP; i++ {
		<-point_done
	}
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

func read_matrix(nrows, ncols int) {
	for i := 0; i < nrows; i++ {
		for j := 0; j < ncols; j++ {
			matrix[i*ncols+j] = byte(read_integer())
		}
	}
}

func read_mask(nrows, ncols int) {
	for i := 0; i < nrows; i++ {
		for j := 0; j < ncols; j++ {
			mask[i][j] = (read_integer() == 1)
		}
	}
}

func FillOnTheFly(nrows, ncols, thresh int) {
  NP := runtime.GOMAXPROCS(0)
  sem := make(chan bool, NP);
  
  var chunkSize int = (nrows + NP - 1) / NP
  var threshInverse = 100 / thresh
  
  // fmt.Printf("%d\n", threshInverse)
  
  for n := 0; n < NP; n++ {
    go func(n int){
    
      var ix int
      var start int = n * chunkSize
      var end int = min( start + chunkSize, nrows)
      
      for i := start; i < end; i++ {
        for j := 0; j < ncols; j++ {
        
          if ((i*ncols+j) % threshInverse ) == 0 {
            mask[i][j] = true 
          }else{
            mask[i][j] = false
          }
          
          ix = i*ncols+j
          matrix[ix] = byte(ix % 100)
        }
      }
      
      sem <- true
    }(n)
  }
  
  // behaves like a barrier
  for n := 0; n < NP; n++ { <-sem }
}

func main() {
	var nrows, ncols, nelts, thresh int
  thresh=0

  flag.Parse()

	nrows = int(read_integer())
	ncols = int(read_integer())

	m := NewByteMatrix(nrows, ncols)

	matrix = m.array

	mask = make([][]bool, nrows)
	for i := range mask {
		mask[i] = make([]bool, ncols)
	}
  
	if !*is_bench {
		read_matrix(nrows, ncols)
		read_mask(nrows, ncols)
	}

	nelts = int(read_integer())
  
  if *is_bench {
    thresh = int(read_integer())
    FillOnTheFly( nrows, ncols, thresh )
  }
  
	points = make([]int, nelts)
  
  var start, stop C.struct_timespec
  var accum C.double
  
  if( C.clock_gettime( C.CLOCK_MONOTONIC_RAW, &start) == -1 ) {
    C.perror( C.CString("clock gettime error 1") );
    return
  }

	Winnow(m, nrows, ncols, nelts)
  
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
  fmt.Fprintf(w, "Go    ,Winnow ,%5d,%5d,%3d,%5d,%2d,%.9f,isBench:%t\n", nrows, ncols, thresh, nelts, runtime.GOMAXPROCS(0), accum, *is_bench )
  
  w.Flush()
  file.Close()

	if !*is_bench {
		fmt.Printf("%d\n", nelts)
		for i := 0; i < nelts; i++ {
			fmt.Printf("%d %d\n", points[i]/ncols, points[i]%ncols)
		}
		fmt.Printf("\n")
	}
}
