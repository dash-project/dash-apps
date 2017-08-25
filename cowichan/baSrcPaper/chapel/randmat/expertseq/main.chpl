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
config const nrows = read(uint (32)),
             ncols = read(uint (32)),
             s     = read(uint (32));

const RowSpace = {1..nrows};
const ColSpace = {1..ncols};
var matrix: [1..nrows, 1..ncols] int (32);


proc randmat() {
  const LCG_A: uint(32) = 1664525,
        LCG_C: uint(32) = 1013904223;
  for i in RowSpace {
    var seed: uint(32) = s + i - 1;
    for j in ColSpace {
      seed = LCG_A * seed + LCG_C;
      matrix[i, j] = (seed % 100) : int (32);
    }
  }
}

proc main() {

  randmat();

  if (!is_bench) {
    writeln(nrows, " ", ncols);

    for i in 1..nrows {
      for j in 1..ncols {
        write(matrix[i, j], " ");
      }
      writeln();
    }
    writeln();
  }
}
