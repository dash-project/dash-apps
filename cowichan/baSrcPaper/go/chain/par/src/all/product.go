/*
 * product: matrix-vector product
 *
 * input:
 *   nelts: the number of elements
 *   Outer_matrix: the real matrix
 *   Outer_vector: the real vector
 *
 * output:
 *   Product_result: a real vector, whose values are the result of the product
 */
package all

var Product_result []Double;

func fill_result_impl(begin, end, ncols int, done chan bool) {
  if (begin + 1 == end) {
    var sum Double = 0;
    for j := 0; j < ncols; j++ {
      sum += Outer_matrix[begin][j] * Outer_vector[j];
    }
    Product_result[begin] = sum;
    done <- true
  } else {
    middle := begin + (end - begin) / 2;
    go fill_result_impl(begin, middle, ncols, done);
    fill_result_impl(middle, end, ncols, done);
  }
}

func fill_result(nrows, ncols int) {
  done := make(chan bool);
  // parallel for on rows
  go fill_result_impl(0, nrows, ncols, done);
  for i := 0; i < nrows; i++ {
    <-done;
  }
}

func Product(nelts int) {
  Product_result = make ([]Double, nelts)
  fill_result(nelts, nelts);
}
