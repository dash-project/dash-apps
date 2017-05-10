/* chain: chain all problems
 * 
 * input:
 *   nelts: the number of elements
 *   randmat_seed: random number generator seed
 *   thresh_percent: percentage of cells to retain
 *   winnow_nelts: the number of points to select
 *
 * output:
 *   result: a real vector, whose values are the result of the final product
 */

use Config, Randmat, Thresh, Winnow, Outer, Product;
config const is_bench = false;
proc main() {
  randmat(nelts, nelts, randmat_seed);
  thresh(nelts, nelts, thresh_percent);
  winnow(nelts, nelts, winnow_nelts);
  outer(winnow_nelts);
  product(winnow_nelts);

  if (!is_bench) {
    writeln(winnow_nelts);
    for i in 1..winnow_nelts do {
      writef("%.4dr ", result[i]);
    }
    writeln();
  }
}
