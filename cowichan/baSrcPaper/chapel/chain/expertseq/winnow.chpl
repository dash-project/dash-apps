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
  var count = 0;
  var values: [0..nrows*ncols] (int, (int, int));
  const RowSpace = {1..nrows};
  const ColSpace = {1..ncols};

  for i in RowSpace {
    for j in ColSpace {
      if (mask[i, j]) {
        values[count] = (matrix[i, j], (i, j));
        count += 1;
      }
    }
  }

  QuickSort(values[0..count]);

  var chunk: int = count / nelts;

  forall i in 1..nelts do {
    var ind: int;
    ind = (i - 1) * chunk + 1;
    (_, points[i]) = values[ind];
  }
}
}
