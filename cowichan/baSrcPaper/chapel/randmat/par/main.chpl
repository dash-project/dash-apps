/* randmat: random number generation
 * 
 * input:
 *   nrows, ncols: the number of rows and columns
 *   s: the seed
 *
 * output:
 *   matrix: an nrows x ncols integer matrix
 */

use Random;

config const is_bench = false;


proc randmat(nrows: int, ncols: int, s: int) {
  var matrix: [1..nrows, 1..ncols]int;
  const LCG_A: uint(32) = 1664525;
  const LCG_C: uint(32) = 1013904223;
  forall i in 1..nrows do {
    var seed: uint(32) = (s + i - 1):uint(32);
    for j in 1..ncols do {
      seed = LCG_A * seed + LCG_C;
      matrix[i, j] = abs(seed) % 100;
    }
  }
  return matrix;
}

proc main() {
  var nrows: int;
  var ncols: int;
  var s: int;

  read(nrows, ncols, s);
  var matrix: [1..nrows, 1..ncols]int;

  matrix = randmat(nrows, ncols, s);

  if (!is_bench) {
    writeln(nrows, " ", ncols);

    for i in 1..nrows do {
      for j in 1..ncols do {
        write(matrix[i, j], " ");
      }
      writeln();
    }
    writeln();
  }
}
