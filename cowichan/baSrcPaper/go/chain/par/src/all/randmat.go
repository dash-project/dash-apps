/*
 * randmat: random number generation
 *
 * input:
 *   nrows, ncols: the number of rows and columns
 *   s: the seed
 *
 * output:
 *   Randmat_matrix: a nrows x ncols integer matrix
 *
 */
package all

var Randmat_matrix [][]byte;

func randvec(row, n, seed int, done chan bool) {
	LCG_A := uint32(1664525)
	LCG_C := uint32(1013904223)
  s := uint32(seed)
  for j := 0; j < n; j++ {
	  s = LCG_A*s + LCG_C
		Randmat_matrix[row][j] = byte(s % 100)
	}
	done <- true
}

func parallel_for(begin, end, ncols, s int, done chan bool) {
  if (begin + 1 == end) {
    randvec(begin, ncols, s + begin, done);
  } else {
    middle := begin + (end - begin) / 2;
    go parallel_for(begin, middle, ncols, s, done);
    parallel_for(middle, end, ncols, s, done);
  }
}

func Randmat(nrows, ncols, s int) {
  Randmat_matrix = make ([][]byte, nrows)
  for i := range Randmat_matrix {
    Randmat_matrix [i] = make ([]byte, ncols)
  }
  
  done := make(chan bool);
  // parallel for on rows
  go parallel_for(0, nrows, ncols, s, done);
  for i := 0; i < nrows; i++ {
    <-done;
  }
}
