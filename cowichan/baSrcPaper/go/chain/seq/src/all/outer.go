/*
 * outer: outer product
 *
 * input:
 *   Outer_vector: a vector of (x, y) points
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

func distance(a, b WinnowPoint) Double {
  return Double(math.Sqrt(float64(sqr(Double(a.i - b.i)) +
        sqr(Double(a.j - b.j)))));
}

func Outer(nelts int) {
  Outer_matrix = make ([][]Double, nelts)
  for i := range Outer_matrix {
    Outer_matrix [i] = make ([]Double, nelts)
  }

  Outer_vector = make ([]Double, nelts)

  for i := 0; i < nelts; i++ {
    var nmax Double = -1;
    for j := 0; j < nelts; j++ {
      if (i != j) {
        Outer_matrix[i][j] = distance(Winnow_points[i], Winnow_points[j]);
        nmax = max(nmax, Outer_matrix[i][j]);
      }
    }
    Outer_matrix[i][i] = Double(nelts) * nmax;
    Outer_vector[i] = distance(WinnowPoint{i: 0, j: 0}, Winnow_points[i]);
  }
}
