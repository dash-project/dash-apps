/** 
 * Measures the latency between each pair
 * of MPI ranks.
 *  
 * To analyze the data see the R scripts
 * provided in the rscript folder
 *
 * author(s): Felix Moessbauer, LMU Munich */

#include <libdash.h>
#include "all-pairs.h"

int main(int argc, char ** argv)
{
  // Gather Parameters
  //
  // Print Setup
  //
  dash::init(&argc, &argv);

  int repeats = 1;

  AllPairs aptest(repeats);
  AllPairsKernel defkern;

  // Run Kernel
  aptest.runKernel(defkern);

  dash::finalize();
}
