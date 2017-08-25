/*
 * outer: outer product
 *
 * input:
 *   vector: a vector of (x, y) points
 *   nelts: the number of points
 *
 * output:
 *   Outer_matrix: a real matrix, whose values are filled with inter-point
 *     distances
 *   Outer_vector: a real vector, whose values are filled with origin-to-point
 *     distances
 */

package all

import (
  "math"
)

type Double float64;

var Outer_matrix [][]Double;
var Outer_vector []Double;

func max(a, b Double) Double {
  if a > b {
    return a;
  }
  return b;
}

func sqr(x Double) Double {
  return x * x;
}

func distance(a, b Point) Double {
  return Double(math.Sqrt(float64(sqr(Double(a.i - b.i)) +
        sqr(Double(a.j - b.j)))));
}

func fill_matrix_impl(begin, end, ncols int, done chan bool) {
  if (begin + 1 == end) {
    var nmax Double = -1;
    for j := 0; j < ncols; j++ {
      if (begin != j) {
        Outer_matrix[begin][j] = distance(Winnow_points[begin], Winnow_points[j]);
        nmax = max(nmax, Outer_matrix[begin][j]);
      }
    }
    Outer_matrix[begin][begin] = Double(ncols) * nmax;
    Outer_vector[begin] = distance(Point{i: 0, j: 0}, Winnow_points[begin]);
    done <- true
  } else {
    middle := begin + (end - begin) / 2;
    go fill_matrix_impl(begin, middle, ncols, done);
    fill_matrix_impl(middle, end, ncols, done);
  }
}

func fill_matrix(nrows, ncols int) {
  done := make(chan bool);
  // parallel for on rows
  go fill_matrix_impl(0, nrows, ncols, done);
  for i := 0; i < nrows; i++ {
    <-done;
  }
}

func Outer(nelts int) {
  Outer_matrix = make ([][]Double, nelts)
  for i := range Outer_matrix {
    Outer_matrix [i] = make ([]Double, nelts)
  }

  Outer_vector = make ([]Double, nelts)

  fill_matrix(nelts, nelts);
}
