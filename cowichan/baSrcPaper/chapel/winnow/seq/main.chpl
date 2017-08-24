/* winnow: weighted point selection
 * 
 * input:
 *   matrix: an integer matrix, whose values are used as masses
 *   mask: a boolean matrix showing which points are eligible for
 *     consideration
 *   nrows, ncols: the number of rows and columns
 *   nelts: the number of points to select
 *
 * output:
 *   points: a vector of (x, y) points
 */

config const is_bench = false;
config const nrows = read(int),
             ncols = read(int);
 
var matrix: [1..nrows, 1..ncols]int;
var mask: [1..nrows, 1..ncols]bool;
var count_per_line: [1..nrows + 1]int;
var values: [0..nrows*ncols] (int, (int, int)); // (value, i, j))

proc winnow(nrows: int, ncols: int, nelts: int) {
  var n: int = 0;
  var points: [1..nelts] (int, int);

  for i in 1..nrows do {
    for j in 1..ncols do {
      if (is_bench) {
        mask[i, j] = (((i - 1) * (j - 1)) % (ncols + 1)) == 1;
      }
      if (mask[i, j]) {
        n += 1;
      }
    }
  }

  var count: int = 1;
  for i in 1..nrows do {
    for j in 1..ncols do {
      if (mask[i, j]) {
        values[count] = (matrix[i, j], (i, j));
        count += 1;
      }
    }
  }

  QuickSort(values[0..n]);

  var chunk: int = n / nelts;

  for i in 1..nelts do {
    var ind: int;
    ind = (i - 1) * chunk + 1;
    (_, points[i]) = values[ind];
  }
  return points;
}

proc read_matrix(nrows, ncols: int) {
  for i in 1..nrows do {
    for j in 1..ncols do {
      read(matrix[i, j]);
    }
  }
}

proc read_mask(nrows, ncols: int) {
  for i in 1..nrows do {
    for j in 1..ncols do {
      var v: int;
      read(v);
      mask[i, j] = v == 1;
    }
  }
}

proc main() {
  var nelts: int;

  if (!is_bench) {
    read_matrix(nrows, ncols);
    read_mask(nrows, ncols);
  }

  read(nelts);
  var points: [1..nelts] (int, int);

  points = winnow(nrows, ncols, nelts);

  if (!is_bench) {
    writeln(nelts);

    for i in 1..nelts do {
      writeln(points[i]);
    }

    writeln();
  }
}
