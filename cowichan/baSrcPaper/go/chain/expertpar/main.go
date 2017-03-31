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
	"math"
	"runtime"
	"sort"
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

func WrapBytes(r, c int, bytes []byte) *ByteMatrix {
	return &ByteMatrix{r, c, bytes}
}

func (m *ByteMatrix) Bytes() []byte {
	return m.array[0 : m.Rows*m.Cols]
}

const (
	LCG_A = 1664525
	LCG_C = 1013904223
)

var (
	is_bench   = flag.Bool("is_bench", false, "")
)

func Randmat(nelts int, s uint32) *ByteMatrix {
	matrix := NewByteMatrix(nelts, nelts)

	work := make(chan int, 256)

	go func() {
		for i := 0; i < nelts; i++ {
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

func Thresh(m *ByteMatrix, nelts, percent int) (mask []bool) {
	NP := runtime.GOMAXPROCS(0)

	hist_work := make(chan int)
	hist_parts := make(chan []int)
	go func() {
		for i := 0; i < nelts; i++ {
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

	count := (nelts * nelts * percent) / 100
	prefixsum := 0
	threshold := 99

	for ; threshold > 0; threshold-- {
		prefixsum += hist[threshold]
		if prefixsum > count {
			break
		}
	}

	mask_work := make(chan int)

	go func() {
		for i := 0; i < nelts; i++ {
			mask_work <- i
		}
		close(mask_work)
	}()

	mask = make([]bool, nelts*nelts)
	mask_done := make(chan bool)
	for i := 0; i < NP; i++ {
		go func() {
			for i := range mask_work {
				row := m.Row(i)
				for j := range row {
					mask[i*nelts+j] = row[j] >= byte(threshold)
				}
			}
			mask_done <- true
		}()
	}

	for i := 0; i < NP; i++ {
		<-mask_done
	}

	return
}

// Winnow structure and sorting helpers
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

type Point struct {
	x, y int
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

func Winnow(m *ByteMatrix, mask []bool, nelts, winnow_nelts int) (points []Point) {
	NP := runtime.GOMAXPROCS(0)
	var values WinnowPoints
	values.m = m

	values_work := make(chan int, 1024)
	values_done := make(chan WinnowPoints, NP)
	values_done <- WinnowPoints{m, make([]int, 0)}

	go func() {
		for i := 0; i < nelts; i++ {
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
				for j := 0; j < nelts; j++ {
					idx := i*nelts + j
					if mask[idx] {
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

	chunk := values.Len() / winnow_nelts

	points = make([]Point, winnow_nelts)
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
				v := values.e[i*chunk]
				points[i] = Point{v / nelts, v % nelts}
			}
			point_done <- true
		}()
	}

	for i := 0; i < NP; i++ {
		<-point_done
	}
	return
}

func Sqr(x float64) float64 {
	return x * x
}

func Distance(ax, ay, bx, by int) float64 {
	return math.Sqrt(float64(Sqr(float64(ax-bx)) + Sqr(float64(ay-by))))
}

func Outer(wp []Point, nelts int) (m [][]float64, vec []float64) {
	m = make([][]float64, nelts)
	vec = make([]float64, nelts)

	NP := runtime.GOMAXPROCS(0)
	work := make(chan int, 256)
	done := make(chan bool)

	go func() {
		for i := range wp {
			work <- i
		}
		close(work)
	}()

	for i := 0; i < NP; i++ {
		go func() {
			for i := range work {
				m[i] = make([]float64, nelts)
				v := wp[i]
				nmax := float64(0)
				for j, w := range wp {
					if i != j {
						d := Distance(v.x, v.y, w.x, w.y)
						if d > nmax {
							nmax = d
						}
						m[i][j] = d
					}
				}
				m[i][i] = float64(nelts) * nmax
				vec[i] = Distance(0, 0, v.x, v.y)
			}
			done <- true
		}()
	}

	for i := 0; i < NP; i++ {
		<-done
	}
	return
}

func Product(m [][]float64, vec []float64, nelts int) (result []float64) {
	result = make([]float64, nelts)
	NP := runtime.GOMAXPROCS(0)
	work := make(chan int)
	done := make(chan bool)

	go func() {
		for i := 0; i < nelts; i++ {
			work <- i
		}
		close(work)
	}()

	for i := 0; i < NP; i++ {
		go func() {
			for i := range work {
				sum := 0.0
				for j := 0; j < nelts; j++ {
					sum += m[i][j] * vec[j]
				}
				result[i] = sum
			}
			done <- true
		}()
	}

	for i := 0; i < NP; i++ {
		<-done
	}
	return

}

func main() {
	flag.Parse()
  var nelts, thresh_percent, seed, winnow_nelts int

	fmt.Scan(&nelts)
	fmt.Scan(&seed)
	fmt.Scan(&thresh_percent)
	fmt.Scan(&winnow_nelts)

	rand_matrix := Randmat(nelts, uint32(seed))
	mask := Thresh(rand_matrix, nelts, thresh_percent)
	win_pts := Winnow(rand_matrix, mask, nelts, winnow_nelts)
	out_matrix, out_vec := Outer(win_pts, winnow_nelts)
	result := Product(out_matrix, out_vec, winnow_nelts)

	if !*is_bench {
		for i := 0; i < winnow_nelts; i++ {
			fmt.Printf("%.3f ", result[i])
		}
		fmt.Printf("\n")
	}
}
