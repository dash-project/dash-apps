/*
 * chain: chain all problems
 *
 * input:
 *   nelts: the number of elements
 *   randmat_seed: random number generator seed
 *   thresh_percent: percentage of cells to retain
 *   winnow_nelts: the number of points to select
 *
 * output:
 *   result: a real vector, whose values are the result of the final product
 *
 */

package main

import (
	"all"
	"flag"
	"fmt"
)

func read_integer() int {
	var value int
	for true {
		var read, _ = fmt.Scanf("%d", &value)
		if read == 1 {
			break
		}
	}
	return value
}

func main() {
	flag.Parse()

	nelts := read_integer()
	randmat_seed := read_integer()
	thresh_percent := read_integer()
	winnow_nelts := read_integer()

	all.Randmat(nelts, nelts, randmat_seed)
	all.Thresh(nelts, nelts, thresh_percent)
	all.Winnow(nelts, nelts, winnow_nelts)
	all.Outer(winnow_nelts)
	all.Product(winnow_nelts)

	if !*all.Is_bench {
		fmt.Printf("%d\n", winnow_nelts)
		for i := 0; i < winnow_nelts; i++ {
			fmt.Printf("%g ", all.Product_result[i])
		}
		fmt.Printf("\n")
	}
}
