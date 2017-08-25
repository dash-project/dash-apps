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

var histogram: [0..99] int;

proc thresh(nrows: int, ncols: int, percent: int) {
  var nmax = max reduce matrix;

  for e in matrix {
    histogram[e] += 1;
  }

  var count: int = (nrows * ncols * percent) / 100;

  var prefixsum: int = 0;
  var threshold: int = nmax;

  for i in histSpace {
    if (prefixsum > count) then break;
    prefixsum += histogram[99 - i];
    threshold = 99 - i;
  }

    mask = matrix >= threshold;
}

}
