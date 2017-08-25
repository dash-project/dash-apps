/* randmat: random number generation
 * 
 * input:
 *   nrows, ncols: the number of rows and columns
 *   s: the seed
 *
 * output:
 *   matrix: an nrows x ncols integer matrix
 */

module Randmat {

use Config;

proc randmat(nrows: int, ncols: int, s: int) {
  const LCG_A: uint(32) = 1664525;
  const LCG_C: uint(32) = 1013904223;

  forall i in 1..nrows do {
    var seed: uint(32) = (s + i - 1): uint(32);
    for j in 1..ncols do {
      seed = LCG_A * seed + LCG_C;
      matrix[i, j] = abs(seed) % 100;
    }
  }
}

}
