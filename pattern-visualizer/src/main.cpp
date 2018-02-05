/**
 * \example ex.06.pattern-block-visualizer/main.cpp
 * Example demonstrating the instantiation of 
 * different patterns and their visualization.
 */ 

#include <main.h>
#include <input_param.h>
#include <pattern_creation.h>
#include <output_pattern.h>
#include <output_info.h>

using namespace dash;

int main(int argc, char* argv[])
{
  dash::init(&argc, &argv);

  auto params = parse_args(argc, argv);
  print_params(params);

  if (dash::myid() == 0) {
    try {
      dash::SizeSpec<2, extent_t> sizespec(params.size_y,  params.size_x);
      dash::TeamSpec<2, index_t>  teamspec(params.units_y, params.units_x);

      if(params.balance_extents) {
        teamspec.balance_extents();
      }
      if (params.tile_x < 0 && params.tile_y < 0) {
        auto max_team_extent = std::max(teamspec.extent(0),
                                        teamspec.extent(1));
        params.tile_y = sizespec.extent(0) / max_team_extent;
        params.tile_x = sizespec.extent(1) / max_team_extent;
      }

      if (params.type == "summa") {
        auto pattern = make_summa_pattern(params, sizespec, teamspec);
        print_example(pattern, params);
      } else if (params.type == "block") {
        auto pattern = make_block_pattern(params, sizespec, teamspec);
        print_example(pattern, params);
      } else if (params.type == "tile") {
        auto pattern = make_tile_pattern(params, sizespec, teamspec);
        print_example(pattern, params);
      } else if (params.type == "shift") {
        auto pattern = make_shift_tile_pattern(params, sizespec, teamspec);
        print_example(pattern, params);
      } else if (params.type == "seq") {
        auto pattern = make_seq_tile_pattern(params, sizespec, teamspec);
        print_example(pattern, params);
      } else {
        print_usage(argv);
        exit(EXIT_FAILURE);
      }

    } catch (std::exception & excep) {
      std::cerr << excep.what() << std::endl;
      if (params.cout) {
        std::cout << "{\"success\": false, \"error\": \"" << excep.what() << "\"}" << std::endl;
      }
    }
  }

  dash::finalize();

  return EXIT_SUCCESS;
}
