#include <boost/program_options.hpp>
#include <libdash.h>
#include <string>
#include <vector>

namespace po = boost::program_options;

/**
* Parses the program options and passes it to the app.
* If everything is correct, valid_opts is set to true
*/
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
    ("ptests", po::value<int>()->default_value(0),
              "number of simultaneously tested pairs. Zero if no limit")
    ("kernels", po::value<std::vector<std::string>>()->multitoken(),
     "kernels to run [def mpi_rma_get mpi_rma_put mpi_sync mpi_async dash_get]")
    ("make_symmetric", po::value<bool>()->default_value(false), "test only upper half plane")
    ("verbose", po::value<int>()->default_value(0), "logging level (0-3), where 0 denotes no logging");

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    try {
      po::notify(vm);
    } catch (po::error& e){
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      if(dash::myid() == 0){
        std::cerr << desc << std::endl;
      }
      valid_opts = false;
    }
  } catch(std::exception& e) {
    std::cerr << "Unhandled Exception reached the top of main: " 
              << e.what() << ", application will now exit" << std::endl; 
    valid_opts = false;
  }

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
