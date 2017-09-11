/* product: matrix-vector product
 * 
 * input:
 *   nelts: the number of elements
 *   matrix: a real matrix
 *   vector: a real vector
 *
 * output:
 *    result: a real vector, whose values are the result of the product
 */

module Product {

use Config;

proc product(nelts: int) {
  forall i in 1..nelts do {
    var sum: real = 0;
    for j in 1..nelts do {
      sum += dists[i, j] * vector[j];
    }
    result[i] = sum;
  }
}

}
