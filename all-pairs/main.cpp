/** 
 * Measures the latency between each pair
 * of MPI ranks.
 *  
 * To analyze the data see the R scripts
 * provided in the rscript folder
 *
 * author(s): Felix Moessbauer, LMU Munich */

#include <libdash.h>
#include <boost/program_options.hpp>
#include "all-pairs.h"
#include "kernel/all-pairs-kernel.h"
#include "kernel/rma-get-kernel.h"
#include "kernel/rma-put-kernel.h"

int  g_repeats, g_ireps;
bool g_make_symmetric;
std::list<std::string> kernels;

void print_help(){
  std::cout << "All Pairs by Felix Moessbauer" << std::endl
            << "Runtime Options:"
            << "--repeats=<number of measurment per pair>" << std::endl
            << "--ireps=<number of repeats per measurement>" << std::endl
            << "--kernels=<komma separated list of kernels>" << std::endl
            << "          possible values are: "
            << "RMA_GET, RMA_PUT, MPI_SYNC, MPI_ASYNC, DASH" << std::endl
            << "--make_symmetric=<0|1> test only upper half plane" <<std::endl;
}

bool parse_options(int & argc, char ** & argv){
  if(dash::myid() == 0){
    print_help();
  }
  return false;
}

int main(int argc, char ** argv)
{
  // Gather Parameters
  //
  // Print Setup
  //
  dash::init(&argc, &argv);

  bool valid_opts = parse_options(argc, argv);

  if(valid_opts){
  AllPairs aptest(g_repeats);

  // Run Kernel
  {
    AllPairsKernel defkern(g_ireps);
    aptest.runKernel(defkern);
  }
  {
    RMAGetKernel rma_get(g_ireps);
    aptest.runKernel(rma_get);
  }
  {
    RMAPutKernel rma_put(g_ireps);
    aptest.runKernel(rma_put);
  }
  }


  dash::finalize();
}
