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
#include "kernel/all-pairs-kernel.h"
#include "kernel/rma-get-kernel.h"
#include "kernel/rma-put-kernel.h"

int main(int argc, char ** argv)
{
  // Gather Parameters
  //
  // Print Setup
  //
  dash::init(&argc, &argv);

  int repeats = 1;

  AllPairs aptest(repeats);

  // Run Kernel
  {
    AllPairsKernel defkern;
    aptest.runKernel(defkern);
  }
  {
    RMAGetKernel rma_get;
    aptest.runKernel(rma_get);
  }
  {
    RMAPutKernel rma_put;
    aptest.runKernel(rma_put);
  }



  dash::finalize();
}
