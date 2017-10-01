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
var histogram[][100]int;

func max_int(a, b int) int {
  if a > b {
    return a;
  }
  return b;
}

func reduce_max_impl(begin, end, ncols int, done chan int) {
  if (begin + 1 == end) {
    res := 0
    for j := 0; j < ncols; j++ {
      res = max_int(res, int(Randmat_matrix[begin][j]));
    }
    done <- res;
  } else {
    middle := begin + (end - begin) / 2;
    go reduce_max_impl(begin, middle, ncols, done);
    reduce_max_impl(middle, end, ncols, done);
  }
}

func reduce_max(nrows, ncols int) int {
  done := make(chan int);
  // parallel for on rows
  go reduce_max_impl(0, nrows, ncols, done);
  res := 0
  for i := 0; i < nrows; i++ {
    res = max_int(res, <-done)
  }
  return res
}

func fill_histogram_impl(begin, end, ncols int, done chan bool) {
  if (begin + 1 == end) {
    for j := 0; j < ncols; j++ {
      histogram[begin][Randmat_matrix[begin][j]]++
    }
    done <- true;
  } else {
    middle := begin + (end - begin) / 2;
    go fill_histogram_impl(begin, middle, ncols, done);
    fill_histogram_impl(middle, end, ncols, done);
  }
}

func fill_histogram(nrows, ncols int) {
  done := make(chan bool);
  // parallel for on rows
  go fill_histogram_impl(0, nrows, ncols, done);
  for i := 0; i < nrows; i++ {
    <-done
  }
}

func merge_histogram_impl(begin, end, n int, done chan bool) {
  if (begin + 1 == end) {
    for i := 1; i < n; i++ {
      histogram[0][begin] += histogram[i][begin]
    }
    done <- true;
  } else {
    middle := begin + (end - begin) / 2;
    go merge_histogram_impl(begin, middle, n, done);
    merge_histogram_impl(middle, end, n, done);
  }
}

func merge_histogram(end, nrows int) {
  done := make(chan bool);
  // parallel for on [0, end)
  go merge_histogram_impl(0, end, nrows, done);
  for i := 0; i < end; i++ {
    <-done
  }
}

func fill_mask_impl(begin, end, ncols, threshold int, done chan bool) {
  if (begin + 1 == end) {
    for j := 0; j < ncols; j++ {
      if (int(Randmat_matrix[begin][j]) >= threshold) {
        Thresh_mask[begin][j] = 1
      }
    }
    done <- true;
  } else {
    middle := begin + (end - begin) / 2;
    go fill_mask_impl(begin, middle, ncols, threshold, done);
    fill_mask_impl(middle, end, ncols, threshold, done);
  }
}

func fill_mask(nrows, ncols, threshold int) {
  done := make(chan bool);
  // parallel for on rows
  go fill_mask_impl(0, nrows, ncols, threshold, done);
  for i := 0; i < nrows; i++ {
    <-done
  }
}

func Thresh(nrows, ncols int, percent int) {
  Thresh_mask = make ([][]byte, nrows)
  for i := range Thresh_mask {
    Thresh_mask [i] = make ([]byte, ncols)
  }
  histogram = make ([][100]int, nrows)

  nmax := reduce_max(nrows, ncols);

  fill_histogram(nrows, ncols);
  merge_histogram(nmax + 1, nrows);

  var count int = (nrows * ncols * percent) / 100;
  var prefixsum int = 0;
  var threshold int = nmax;

  for i := nmax; i >= 0 && prefixsum <= count; i-- {
    prefixsum += histogram[0][i];
    threshold = i;
  }

  fill_mask(nrows, ncols, threshold)
}
