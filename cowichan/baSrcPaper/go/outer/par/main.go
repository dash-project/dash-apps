/*
 * outer: outer product
 *
 * input:
 *   vector: a vector of (x, y) points
 *   nelts: the number of points
 *
 * output:
 *   matrix: a real matrix, whose values are filled with inter-point
 *     distances
 *   vector: a real vector, whose values are filled with origin-to-point
 *     distances
 */

package main

import (
	"flag"
	"fmt"
	"math"
)

var is_bench = flag.Bool("is_bench", false, "")

type Point struct {
	i, j int
}

var points []Point

type double float64

var matrix [][]double
var vector []double

type Points []Point

func max(a, b double) double {
	if a > b {
		return a
	}
	return b
}

func sqr(x double) double {
	return x * x
}

func distance(a, b Point) double {
	return double(math.Sqrt(float64(sqr(double(a.i-b.i)) +
		sqr(double(a.j-b.j)))))
}

func fill_values_impl(begin, end, ncols int, done chan bool) {
	if begin+1 == end {
		var nmax double = -1
		for j := 0; j < ncols; j++ {
			if begin != j {
				matrix[begin][j] = distance(points[begin], points[j])
				nmax = max(nmax, matrix[begin][j])
			}
		}
		matrix[begin][begin] = double(ncols) * nmax
		vector[begin] = distance(Point{0, 0}, points[begin])
		done <- true
	} else {
		middle := begin + (end-begin)/2
		go fill_values_impl(begin, middle, ncols, done)
		fill_values_impl(middle, end, ncols, done)
	}
}

func fill_values(nrows, ncols int) {
	done := make(chan bool)
	// parallel for on rows
	go fill_values_impl(0, nrows, ncols, done)
	for i := 0; i < nrows; i++ {
		<-done
	}
}

func outer(nelts int) {
	fill_values(nelts, nelts)
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

func read_vector_of_points(nelts int) {
	for i := 0; i < nelts; i++ {
		a := read_integer()
		b := read_integer()
		points[i] = Point{a, b}
	}
}

func main() {
	var nelts int

	flag.Parse()

	nelts = read_integer()
  matrix = make ([][]double, nelts)
  for i := range matrix {
    matrix [i] = make ([]double, nelts)
  }

  vector = make ([]double, nelts)
  points = make ([]Point, nelts)


	if !*is_bench {
		read_vector_of_points(nelts)
	}

	outer(nelts)

	if !*is_bench {
		fmt.Printf("%d %d\n", nelts, nelts)
		for i := 0; i < nelts; i++ {
			for j := 0; j < nelts; j++ {
				fmt.Printf("%g ", matrix[i][j])
			}
			fmt.Printf("\n")
		}
		fmt.Printf("\n")

		fmt.Printf("%d\n", nelts)
		for i := 0; i < nelts; i++ {
			fmt.Printf("%g ", vector[i])
		}
		fmt.Printf("\n")
	}
}
