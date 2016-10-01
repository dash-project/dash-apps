/**
 * Measures the latency between each pair
 * of MPI ranks.
 *
 * To analyze the data see the R scripts
 * provided in the rscript folder
 *
 * author(s): Felix Moessbauer, LMU Munich */


#include <libdash.h>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include "logger.h"
#include "program_options.h"
#include "all-pairs.h"
#include "kernel/all-pairs-kernel.h"
#include "kernel/rma-get-kernel.h"
#include "kernel/rma-put-kernel.h"
#include "kernel/mpi-sync-kernel.h"
#include "kernel/mpi-async-kernel.h"
#include "kernel/dash-get-kernel.h"

#include <vector>
#include <string>
#include <algorithm>

using kernels_type = std::vector<std::string>;

namespace logging  = boost::log;
namespace logtriv  = logging::trivial;

void setupLogger(int loglevel){
  auto logger = logging::core::get();

  switch(loglevel){
    case 1:
      logger->set_filter(logtriv::severity >= logtriv::info);
      break;
    case 2:
      logger->set_filter(logtriv::severity >= logtriv::debug);
      break;
    case 3:
      logger->set_filter(logtriv::severity >= logtriv::trace);
      break;
    default:
      logger->set_logging_enabled(false);
  }
  LOG_UNIT(info) << "Logging enabled";
}

int main(int argc, char ** argv)
{
    dash::init(&argc, &argv);

    bool valid_opts = true;
    auto opts = setup_program_options(argc, argv, valid_opts);

    if(valid_opts) {
        int  repeats  = opts["repeats"].as<int>();
        int  ireps    = opts["ireps"].as<int>();
        int  ptests   = opts["ptests"].as<int>();
        bool make_sym = opts["make_symmetric"].as<bool>();
        auto kernels  = opts["kernels"].as<kernels_type>();
        int  loglevel = opts["verbose"].as<int>();

        // Sanitize
        if((ptests <= 0) || (ptests > dash::size())){
          ptests = dash::size();
        }

        setupLogger(loglevel);

        AllPairs aptest(repeats, ptests, make_sym);

        for(auto k:kernels) {
            if(k == "def") {
                AllPairsKernel defkern(ireps);
                aptest.runKernel(defkern);
            } else if(k == "mpi_rma_get") {
                RMAGetKernel rma_get(ireps);
                aptest.runKernel(rma_get);
            } else if(k == "mpi_rma_put") {
                RMAPutKernel rma_put(ireps);
                aptest.runKernel(rma_put);
            } else if(k == "mpi_sync") {
                MPISyncKernel mpi_sync(ireps);
                aptest.runKernel(mpi_sync);
            } else if(k == "mpi_async") {
                MPIASyncKernel mpi_async(ireps);
                aptest.runKernel(mpi_async);
            } else if(k == "dash_get") {
                DashGetKernel dash_get(ireps);
                aptest.runKernel(dash_get);
            } else {
                std::cout << "unknown kernel" << std::endl;
            }
        }
    }

    dash::finalize();
}
