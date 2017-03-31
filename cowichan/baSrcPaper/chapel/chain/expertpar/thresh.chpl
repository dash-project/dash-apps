/* thresh: histogram thresholding
 * 
 * input:
 *   matrix: the integer matrix to be thresholded
 *   nrows, ncols: the number of rows and columns
 *   percent: the percentage of cells to retain
 *
 * output:
 *   mask: a boolean matrix filled with true for the cells to be kept
 */

module Thresh {
use Config;

proc thresh(nrows: int, ncols: int, percent: int) {
  var nmax = max reduce matrix;
  var histogram: [1..nrows, 0..99] int;
  const RowSpace = {1..nrows};
  const ColSpace = {1..ncols};

  // it was recommended to use a 1D array of atomic ints,
  // but it was much slower
  forall i in RowSpace {
    for j in ColSpace {
      histogram[i, matrix[i, j]] += 1;
    }
  }

  const RowSpace2 = {2..nrows};

  forall j in 0..(nmax) {
    for i in RowSpace2 {
      histogram[1, j] += histogram[i, j];
    }
  }

  var count: int = (nrows * ncols * percent) / 100;

  var prefixsum: int = 0;
  var threshold: int = nmax;

  for i in histSpace do {
    if (prefixsum > count) then break;
    prefixsum += histogram[1, 99 - i];
    threshold = 99 - i;
  }

  forall i in RowSpace {
    for j in ColSpace {
      mask[i, j] = matrix[i, j] >= threshold;
    }
  }

  // this was recommended, but turned out to be quite a bit slower.
  // mask = matrix >= threshold;
}

}
