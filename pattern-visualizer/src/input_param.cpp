#include <input_param.h>
#include <main.h>

#include <cstdlib>
#include <iostream>
#include <fstream>

void print_basic_params() {
  std::ifstream p("params.json");

  if(p.is_open()) {
    std::cout << p.rdbuf();
  }
}

void print_usage(char **argv)
{
  if (dash::myid() == 0) {
    std::cerr << "Usage: \n"
              << basename(argv[0]) << " "
              << "-h | "
              << "-p | "
              << "[-b] "
              << "[-m] "
              << "[-r]"
              << "\n\n";
    std::cerr << "-h  Print this help and exit\n"
              << "-p  Print the basic structure for pattern parameters and exit\n\n"
              << "Parameters to construct the pattern are read from stdin\n"
              << "-b  Just output the blocks of the pattern, may fail for irregular patterns\n"
              << "-m  Just output the memory layout, coordinates ordered by local offset\n"
              << "-r  Reduced output, incomplete but enough for a first draw"
              << std::endl;
  }
}

/*void print_params(const cli_params & params)
{
  int w_size  = std::log10(std::max(params.size_x,  params.size_y));
  int w_units = std::log10(std::max(params.units_x, params.units_y));
  int w_tile  = std::log10(std::max(params.tile_x,  params.tile_y));
  int w       = std::max({ w_size, w_units, w_tile }) + 1;

  std::cerr << "Parameters:"
            << '\n'
            << "    type (-s):               " << params.type
            << '\n'
            << "    size (-n <rows> <cols>): ( "
            << std::fixed << std::setw(w) << params.size_y << ", "
            << std::fixed << std::setw(w) << params.size_x << " )"
            << '\n'
            << "    team (-u <rows> <cols>): ( "
            << std::fixed << std::setw(w) << params.units_y << ", "
            << std::fixed << std::setw(w) << params.units_x << " )"
            << '\n'
            << "    balance extents (-e): "
            << (params.balance_extents ? "yes" : "no")
            << '\n'
            << "    tile (-t <rows> <cols>): ( "
            << std::fixed << std::setw(w) << params.tile_y << ", "
            << std::fixed << std::setw(w) << params.tile_x << " )"
            << '\n'
            << "    blocked display (-b): "
            << (params.blocked_display ? "yes" : "no")
            << '\n'
            << std::endl;
}*/

cli_params parse_args(int argc, char * argv[])
{
  cli_params params = default_output_params;

  bool resetted = false;
  for (auto i = 1; i < argc; i++) {
    std::string flag = argv[i];
    if (flag == "-h") {
      print_usage(argv);
      exit(EXIT_SUCCESS);
    }
    if (flag == "-p") {
      print_basic_params();
      exit(EXIT_SUCCESS);
    }
    if (flag == "-b" || flag == "-m") {
      if(!resetted) {
        params.output_blocks = false;
        params.output_memlayout = false;
        resetted = true;
      }
    }
    if (flag == "-b") {
      params.output_blocks = true;
    } else if (flag == "-t") {
      params.output_tiles = true;
    } else if (flag == "-m") {
      params.output_memlayout = true;
    } else if (flag == "-r") {
      params.reduced = true;
    } else {
      print_usage(argv);
      exit(EXIT_FAILURE);
    }
  }

  return params;
}

json read_params() {
  json params;

  std::cin >> params;

  return params;
}

const json getParamFromGroup(const json contents, const std::string name) {
  return *std::find_if(contents.begin(), contents.end(),
                       [&](const json &m){
                         return m.at("name") == name;
                       });
}

/*json getParamsFromSelection(json selection, json value) {
  if(value.is_number()) {
    return selection.at("values").at(value).at("params");
  }
  // todo
  // todo quantity of selections in groups
}*/
