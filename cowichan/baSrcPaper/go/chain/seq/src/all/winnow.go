/*
 * winnow: weighted point selection
 *
 * input:
 *   Randmat_matrix: an integer matrix, whose values are used as masses
 *   Thresh_mask: a boolean matrix showing which points are eligible for
 *     consideration
 *   nrows, ncols: the number of rows and columns
 *   nelts: the number of points to select
 *
 * output:
 *   Winnow_points: a vector of (x, y) points
 *
 */
package all

import (
  "sort"
  "flag"
)

var Is_bench = flag.Bool("is_bench", false, "")

type WinnowPoint struct {
  value byte;
  i, j int;
}

var Winnow_points []WinnowPoint;

type WinnowPoints []WinnowPoint;

func (p WinnowPoints) Len() int { return len(p) }
func (p WinnowPoints) Swap(i, j int) { p[i], p[j] = p[j], p[i] }
func (p WinnowPoints) Less(i, j int) bool {
  if p[i].value < p[j].value {
    return true;
  }
  if p[i].value > p[j].value {
    return false;
  }
  if p[i].i < p[j].i {
    return true;
  }
  if p[i].i > p[j].i {
    return false;
  }
  return p[i].j < p[j].j;
}

func Winnow(nrows, ncols, nelts int) {
  Winnow_points = make ([]WinnowPoint, nelts)
  var n = 0;

  for i := 0; i < nrows; i++ {
    for j := 0; j < ncols; j++ {
      if (*Is_bench) {
        if (((i * j) % (ncols + 1)) == 1) {
          Thresh_mask[i][j] = 1;
        }
      }
      if Thresh_mask[i][j] == 1 {
        n++;
      }
    }
  }

  values := make(WinnowPoints, n);

  var count = 0;
  for i := 0; i < nrows; i++ {
    for j := 0; j < ncols; j++ {
      if Thresh_mask[i][j] == 1 {
        values[count] = WinnowPoint{Randmat_matrix[i][j], i, j};
        count++;
      }
    }
  }
  
  sort.Sort(values);

  var total = len(values);
  var chunk int = total / nelts;

  for i := 0; i < nelts; i++ {
    var index = i * chunk;
    Winnow_points[i] = values[index];
  }
}
