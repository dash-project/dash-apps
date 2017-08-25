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

config const is_bench = false;
config const nelts = read(int);
var matrix: [1..nelts, 1..nelts]real;
var vector: [1..nelts]real;
var result: [1..nelts]real;

proc product(nelts: int) {
  forall i in 1..nelts do {
    var sum: real = 0;
    for j in 1..nelts do {
      sum += matrix[i, j] * vector[j];
    }
    result[i] = sum;
  }
}

proc main() {
  if (!is_bench) {
    for i in 1..nelts do {
      for j in 1..nelts do {
        read(matrix[i, j]);
      }
    }

    for i in 1..nelts do {
      read(vector[i]);
    }
  }

  product(nelts);

  if (!is_bench) {
    writeln(nelts);
    for i in 1..nelts do {
      write(result[i] + " ");
    }
    writeln();
  }
}
