/* thresh: histogram thresholding
 * 
 * input:
 *   matrix: the integer matrix to be thresholded
 *   nrows, ncols: the number of rows and columns
 *   percent: the percentage of cells to retain
 *
 * output:
 *   mask: a boolean matrix filled with true for the cells to be kept
 */

config const is_bench = false;
config const nrows = read(int),
             ncols = read(int);

const ProbSpace = {1..nrows, 1..ncols},
      HistSpace = {0..100};
var matrix: [ProbSpace] int; 
var mask: [ProbSpace] int;
var histogram: [HistSpace] int;

proc thresh(nrows: int, ncols: int, percent: int) {
  var nmax = max reduce matrix;

  for (i,j) in ProbSpace {
    histogram[matrix[i,j]] += 1;
  }

  var count: int = (nrows * ncols * percent) / 100;

  var prefixsum: int = 0;
  var threshold: int = nmax;

  for i in 0..100 {
    if (prefixsum > count) then break;
    prefixsum += histogram[100 - i];
    threshold = 100 - i ;
  }

  for (i, j) in ProbSpace {
    mask[i, j] = matrix[i, j] >= threshold;
  }
}

proc main() {
  var percent: int;

  if (!is_bench) {
    for i in 1..nrows {
      for j in 1..ncols {
        read(matrix[i,j]);
      }
    }
  }

  read(percent);

  thresh(nrows, ncols, percent);

  if (!is_bench) {
    writeln(nrows, " ", ncols);

    for i in 1..nrows {
      for j in 1..ncols {
        if (mask[i, j]) {
          write("1", " ");
        } else {
          write("0", " ");
        }
      }
      writeln();
    }
    writeln();
  }
}
