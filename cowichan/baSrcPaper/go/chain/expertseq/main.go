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
	is_bench = flag.Bool("is_bench", false, "")
)

func Randmat(nelts int, s uint32) *ByteMatrix {
	matrix := NewByteMatrix(nelts, nelts)

	for i := 0; i < nelts; i++ {
		var seed = s + uint32(i)
		row := matrix.Row(i)
		for j := range row {
			seed = LCG_A*seed + LCG_C
			row[j] = byte(seed%100) % 100
		}
	}
	return matrix
}

func Thresh(m *ByteMatrix, nelts, percent int) (mask []bool) {
	var hist [100]int
	mask = make([]bool, nelts*nelts)
	for _, v := range m.Bytes() {
		hist[v]++
	}

	count := (nelts * nelts * percent) / 100
	prefixsum := 0
	var threshold int

	for threshold = 99; threshold > 0; threshold-- {
		prefixsum += hist[threshold]
		if prefixsum > count {
			break
		}
	}

	for i := 0; i < nelts; i++ {
		row := m.Row(i)
		for j := range row {
			mask[i*nelts+j] = row[j] >= byte(threshold)
		}
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
	if p.m.array[p.e[i]] != p.m.array[p.e[j]] {
		return p.m.array[p.e[i]] < p.m.array[p.e[j]]
	}

	return p.e[i] < p.e[j]
}

type Point struct {
	x, y int
}

func Winnow(m *ByteMatrix, mask []bool, nelts, winnow_nelts int) (points []Point) {
	var values WinnowPoints
	values.m = m

	for i := 0; i < nelts; i++ {
		for j := 0; j < nelts; j++ {
			idx := i*nelts + j
			if mask[idx] {
				values.e = append(values.e, idx)
			}
		}
	}
	sort.Sort(&values)

	chunk := values.Len() / winnow_nelts

	points = make([]Point, winnow_nelts)
	for i := 0; i < winnow_nelts; i++ {
		v := values.e[i*chunk]
		p := Point{v / nelts, v % nelts}
		points[i] = p
	}

	return
}

func Sqr(x float64) float64 {
	return x * x
}

func Distance(ax, ay, bx, by int) float64 {
	return math.Sqrt(float64(Sqr(float64(ax-bx)) + Sqr(float64(ay-by))))
}

func Outer(wp []Point, nelts int) (m []float64, vec []float64) {
	m = make([]float64, nelts*nelts)
	vec = make([]float64, nelts)
	for i, v := range wp {
		nmax := float64(0)
		for j, w := range wp {
			if i != j {
				d := Distance(v.x, v.y, w.x, w.y)
				if d > nmax {
					nmax = d
				}
				m[i*nelts+j] = d
			}
		}
		m[i*(nelts+1)] = float64(nelts) * nmax
		vec[i] = Distance(0, 0, v.x, v.y)
	}
	return
}

func Product(m, vec []float64, nelts int) (result []float64) {
	result = make([]float64, nelts)
	for i := 0; i < nelts; i++ {
		sum := 0.0
		for j := 0; j < nelts; j++ {
			sum += m[i*nelts+j] * vec[j]
		}
		result[i] = sum
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
