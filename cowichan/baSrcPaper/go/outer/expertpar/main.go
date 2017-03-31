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
	"runtime"
)

var is_bench = flag.Bool("is_bench", false, "")

var points []Point

type Point struct {
	x, y int
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
	work := make(chan int)

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
  points = make ([]Point, nelts)

	if !*is_bench {
		read_vector_of_points(nelts)
	}

	matrix, vector := Outer(points[0:nelts], nelts)

	if !*is_bench {
		fmt.Printf("%d %d\n", nelts, nelts)
    for _, row := range matrix {
      for _, elem := range row {
        fmt.Printf("%g ", elem)
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
