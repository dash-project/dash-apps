/*
 * thresh: histogram thresholding
 *
 * input:
 *   Randmat_matrix: the integer matrix to be thresholded
 *   nrows, ncols: the number of rows and columns
 *   percent: the percentage of cells to retain
 *
 * output:
 *   Thresh_mask: a boolean matrix filled with true for cells that are kept
 *
 */
package all

var Thresh_mask [][]byte;
var histogram [100]int;

func max_int(a, b int) int {
  if a > b {
    return a;
  }
  return b;
}

func Thresh(nrows, ncols int, percent int) {
  Thresh_mask = make ([][]byte, nrows)
  for i := range Thresh_mask {
    Thresh_mask [i] = make ([]byte, ncols)
  }

  var nmax int = 0;
  for i := 0; i < nrows; i++ {
    for j := 0; j < ncols; j++ {
      nmax = max_int(nmax, int(Randmat_matrix[i][j]));
    }
  }

  for i := 0; i < nrows; i++ {
    for j := 0; j < ncols; j++ {
      histogram[Randmat_matrix[i][j]]++;
    }
  }

  var count int = (nrows * ncols * percent) / 100;
  var prefixsum int = 0;
  var threshold int = nmax;

  for i := nmax; i >= 0 && prefixsum <= count; i-- {
    prefixsum += histogram[i];
    threshold = i;
  }

  for i := 0; i < nrows; i++ {
    for j := 0; j < ncols; j++ {
      if int(Randmat_matrix[i][j]) >= threshold {
        Thresh_mask[i][j] = 1;
      }
    }
  }
}
