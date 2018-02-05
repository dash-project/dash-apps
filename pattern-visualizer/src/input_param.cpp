#include <input_param.h>
#include <main.h>

#include <cstdlib>

void print_usage(char **argv)
{
  if (dash::myid() == 0) {
    std::cerr << "Usage: \n"
              << basename(argv[0]) << " "
              << "-h | "
              << "[-s pattern] "
              << "[-n size_spec] "
              << "[-u unit_spec] "
              << "[-t tile_spec] "
              << "[-p] "
              << "\n\n";
    std::cerr << "-s pattern:   [summa|block|tile|seq|shift]\n"
              << "-n size_spec: <size_y>  <size_x>  [ "
              << default_params.size_y << " " << default_params.size_x << " ]\n"
              << "-u unit_spec: <units_y> <units_x> [  "
              << default_params.units_x << "  " << default_params.units_x << " ]\n"
              << "-t tile_spec: <tile_y>  <tile_x>  [ automatically determined ]\n"
              << "-p          : print to stdout instead of stderr\n"
              << "-h          : print help and exit"
              << std::endl;
  }
}

void print_params(const cli_params & params)
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
}

cli_params parse_args(int argc, char * argv[])
{
  cli_params params = default_params;

  for (auto i = 1; i < argc; i += 3) {
    std::string flag = argv[i];
    if (flag == "-h") {
      print_usage(argv);
      exit(EXIT_SUCCESS);
    }
    if (flag == "-s") {
      params.type    = argv[i+1];
      i -= 1;
    } else if (flag == "-n") {
      params.size_y  = static_cast<extent_t>(atoi(argv[i+1]));
      params.size_x  = static_cast<extent_t>(atoi(argv[i+2]));
    } else if (flag == "-u") {
      params.units_y = static_cast<extent_t>(atoi(argv[i+1]));
      params.units_x = static_cast<extent_t>(atoi(argv[i+2]));
    } else if (flag == "-t") {
      params.tile_y  = static_cast<extent_t>(atoi(argv[i+1]));
      params.tile_x  = static_cast<extent_t>(atoi(argv[i+2]));
    } else if (flag == "-p") {
      params.cout    = true;
      i -= 2;
    } else if (flag == "-e") {
      params.balance_extents = true;
      i -= 2;
    } else if (flag == "-b") {
      params.blocked_display = true;
      i -= 2;
    } else {
      print_usage(argv);
      exit(EXIT_FAILURE);
    }
  }

  return params;
}
