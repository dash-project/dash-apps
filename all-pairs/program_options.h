#include <boost/program_options.hpp>
#include <string>
#include <vector>

namespace po = boost::program_options;

po::variables_map setup_program_options(int &argc, char ** &argv, bool &valid_opts){
  if(dash::myid() == 0){
    std::cout << "All Pairs by Felix Moessbauer" << std::endl 
            << "This program measures the latency between all mpi ranks"
            << " using various kernels."
            << std::endl;
  }

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "show help message")
    ("repeats", po::value<int>()->default_value(50), "number of measurements per pair")
    ("ireps", po::value<int>()->default_value(5), "number of repeats per measurement")
    ("kernels", po::value<std::vector<std::string>>()->multitoken(),
     "kernels to run [def mpi_rma_get mpi_rma_put mpi_sync mpi_async dash]")
    ("make_symmetric", po::value<bool>()->default_value(false), "test only upper half plane");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    

  if(vm.count("help") || argc == 1) {
    if(dash::myid() == 0){
      std::cout << desc << std::endl;
    }
    valid_opts = false;
  } else if(vm.count("kernels") == 0){
    if(dash::myid() == 0){
      std::cout << "select at least one kernel" << std::endl;
    }
    valid_opts = false;
  }

  return vm;
}
