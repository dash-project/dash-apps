/* winnow: weighted point selection
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
 */

module Winnow {
use Config;

proc winnow(nrows: int, ncols: int, nelts: int) {
  var n: int = 0;
  var count_per_line: [1..nrows+1] int;
  var values: [0..nrows*ncols] (int, (int, int)); // (value, i, j))
  const RowSpace = {1..nrows};
  const ColSpace = {1..ncols};

  forall i in RowSpace {
    count_per_line[i + 1] = 0;
    for j in ColSpace {
     count_per_line[i + 1] += mask[i, j];
    }
  }

  var total = + scan count_per_line;
  n = total[nrows + 1];

  forall i in RowSpace {
    var count = total[i];
    for j in ColSpace {
      if (mask[i, j]) {
        values[count] = (matrix[i, j], (i, j));
        count += 1;
      }
    }
  }

  QuickSort(values[0..n]);

  var chunk: int = n / nelts;

  forall i in 1..nelts do {
    var ind: int;
    ind = (i - 1) * chunk + 1;
    (_, points[i]) = values[ind];
  }
}
}
