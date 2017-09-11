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

proc product(nelts: int)
{
  for (i, j) in {1..nelts, 1..nelts} {
    result [i] += dists[i, j] * vector[j];
  }
}

}
