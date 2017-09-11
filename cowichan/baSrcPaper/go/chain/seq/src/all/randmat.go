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

func Randmat(nrows, ncols, s int) {
  Randmat_matrix = make ([][]byte, nrows)
  for i := range Randmat_matrix {
    Randmat_matrix [i] = make ([]byte, ncols)
  }
	LCG_A := uint32(1664525)
	LCG_C := uint32(1013904223)
	for i := 0; i < nrows; i++ {
		seed := uint32(s) + uint32(i)
		for j := 0; j < ncols; j++ {
			seed = LCG_A*seed + LCG_C
			Randmat_matrix[i][j] = byte(seed % 100)
		}
	}
}
