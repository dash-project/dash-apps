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

func outer(nelts int) {
	for i := 0; i < nelts; i++ {
		var nmax double = -1
		for j := 0; j < nelts; j++ {
			if i != j {
				matrix[i][j] = distance(points[i], points[j])
				nmax = max(nmax, matrix[i][j])
			}
		}
		matrix[i][i] = double(nelts) * nmax
		vector[i] = distance(Point{0, 0}, points[i])
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
