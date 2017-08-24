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

func Product(nelts int) {
  Product_result = make ([]Double, nelts)
  for i := 0; i < nelts; i++ {
    var sum Double = 0;
    for j := 0; j < nelts; j++ {
      sum += Outer_matrix[i][j] * Outer_vector[j];
    }
    Product_result[i] = sum;
  }
}
