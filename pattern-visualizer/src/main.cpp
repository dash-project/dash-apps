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

  auto output_params = parse_args(argc, argv);

  if (dash::myid() == 0) {
    try {
      const auto params = read_params().at("content")[0];

      if(true){//params.at("value").at("name") == "con") {
        auto con_params = params.at("value").at("params");
        std::string       type_r = getParamFromGroup(con_params,"pattern").at("value").at("name");
        std::string       ndim_str = getParamFromGroup(con_params,"numDim").at("value");
        dash::dim_t       ndim_r = std::stoi(ndim_str);
        dash::MemArrange  arr_r = dash::ROW_MAJOR;
        if(getParamFromGroup(con_params,"arrangement").at("value").at("name") == "col") {
          arr_r = dash::COL_MAJOR;
        }

        if(arr_r == dash::ROW_MAJOR) {
          constexpr dash::MemArrange arr = dash::ROW_MAJOR;
          if(ndim_r == 1) {
            constexpr dash::dim_t ndim = 1;
            if(type_r == "block") {
              auto pattern = make_constructor_pattern<block,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "tile") {
              auto pattern = make_constructor_pattern<tile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "shiftTile") {
              auto pattern = make_constructor_pattern<shiftTile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "seqTile") {
              auto pattern = make_constructor_pattern<seqTile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            /*} else if(type_r == "dynamic") {
              auto pattern = make_constructor_pattern<dynamic,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);*/
            } else if(type_r == "csr") {
              auto pattern = make_constructor_pattern<csr,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            }
          } else
          if(ndim_r == 2) {
            constexpr dash::dim_t ndim = 2;
            if(type_r == "block") {
              auto pattern = make_constructor_pattern<block,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "tile") {
              auto pattern = make_constructor_pattern<tile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "shiftTile") {
              auto pattern = make_constructor_pattern<shiftTile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "seqTile") {
              auto pattern = make_constructor_pattern<seqTile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            /*} else if(type_r == "dynamic") {
              auto pattern = make_constructor_pattern<dynamic,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);*/
            /*} else if(type_r == "csr") {
              auto pattern = make_constructor_pattern<csr,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);*/
            }
          } else if(ndim_r == 3) {
            constexpr dash::dim_t ndim = 3;
            if(type_r == "block") {
              auto pattern = make_constructor_pattern<block,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "tile") {
              auto pattern = make_constructor_pattern<tile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "shiftTile") {
              auto pattern = make_constructor_pattern<shiftTile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "seqTile") {
              auto pattern = make_constructor_pattern<seqTile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            /*} else if(type_r == "dynamic") {
              auto pattern = make_constructor_pattern<dynamic,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);*/
            /*} else if(type_r == "csr") {
              auto pattern = make_constructor_pattern<csr,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);*/
            }
          } else if(ndim_r == 4) {
            constexpr dash::dim_t ndim = 4;
            if(type_r == "block") {
              auto pattern = make_constructor_pattern<block,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "tile") {
              auto pattern = make_constructor_pattern<tile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "shiftTile") {
              auto pattern = make_constructor_pattern<shiftTile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "seqTile") {
              auto pattern = make_constructor_pattern<seqTile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            /*} else if(type_r == "dynamic") {
              auto pattern = make_constructor_pattern<dynamic,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);*/
            /*} else if(type_r == "csr") {
              auto pattern = make_constructor_pattern<csr,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);*/
            }
          } else if(ndim_r == 5) {
            constexpr dash::dim_t ndim = 5;
            if(type_r == "block") {
              auto pattern = make_constructor_pattern<block,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "tile") {
              auto pattern = make_constructor_pattern<tile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "shiftTile") {
              auto pattern = make_constructor_pattern<shiftTile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "seqTile") {
              auto pattern = make_constructor_pattern<seqTile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            /*} else if(type_r == "dynamic") {
              auto pattern = make_constructor_pattern<dynamic,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);*/
            /*} else if(type_r == "csr") {
              auto pattern = make_constructor_pattern<csr,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);*/
            }
          }
        } else if(arr_r == dash::COL_MAJOR) {
          constexpr dash::MemArrange arr = dash::COL_MAJOR;
          if(ndim_r == 1) {
            constexpr dash::dim_t ndim = 1;
            if(type_r == "block") {
              auto pattern = make_constructor_pattern<block,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "tile") {
              auto pattern = make_constructor_pattern<tile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "shiftTile") {
              auto pattern = make_constructor_pattern<shiftTile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "seqTile") {
              auto pattern = make_constructor_pattern<seqTile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            /*} else if(type_r == "dynamic") {
              auto pattern = make_constructor_pattern<dynamic,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);*/
            } else if(type_r == "csr") {
              auto pattern = make_constructor_pattern<csr,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            }
          } else
          if(ndim_r == 2) {
            constexpr dash::dim_t ndim = 2;
            if(type_r == "block") {
              auto pattern = make_constructor_pattern<block,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "tile") {
              auto pattern = make_constructor_pattern<tile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "shiftTile") {
              auto pattern = make_constructor_pattern<shiftTile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "seqTile") {
              auto pattern = make_constructor_pattern<seqTile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            /*} else if(type_r == "dynamic") {
              auto pattern = make_constructor_pattern<dynamic,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);*/
            /*} else if(type_r == "csr") {
              auto pattern = make_constructor_pattern<csr,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);*/
            }
          } else if(ndim_r == 3) {
            constexpr dash::dim_t ndim = 3;
            if(type_r == "block") {
              auto pattern = make_constructor_pattern<block,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "tile") {
              auto pattern = make_constructor_pattern<tile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "shiftTile") {
              auto pattern = make_constructor_pattern<shiftTile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "seqTile") {
              auto pattern = make_constructor_pattern<seqTile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            /*} else if(type_r == "dynamic") {
              auto pattern = make_constructor_pattern<dynamic,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);*/
            /*} else if(type_r == "csr") {
              auto pattern = make_constructor_pattern<csr,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);*/
            }
          } else if(ndim_r == 4) {
            constexpr dash::dim_t ndim = 4;
            if(type_r == "block") {
              auto pattern = make_constructor_pattern<block,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "tile") {
              auto pattern = make_constructor_pattern<tile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "shiftTile") {
              auto pattern = make_constructor_pattern<shiftTile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "seqTile") {
              auto pattern = make_constructor_pattern<seqTile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            /*} else if(type_r == "dynamic") {
              auto pattern = make_constructor_pattern<dynamic,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);*/
            /*} else if(type_r == "csr") {
              auto pattern = make_constructor_pattern<csr,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);*/
            }
          } else if(ndim_r == 5) {
            constexpr dash::dim_t ndim = 5;
            if(type_r == "block") {
              auto pattern = make_constructor_pattern<block,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "tile") {
              auto pattern = make_constructor_pattern<tile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "shiftTile") {
              auto pattern = make_constructor_pattern<shiftTile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            } else if(type_r == "seqTile") {
              auto pattern = make_constructor_pattern<seqTile,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);
            /*} else if(type_r == "dynamic") {
              auto pattern = make_constructor_pattern<dynamic,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);*/
            /*} else if(type_r == "csr") {
              auto pattern = make_constructor_pattern<csr,ndim,arr>(con_params);
              print_example(pattern,output_params,con_params);*/
            }
          }
        }
      }
    } catch (std::exception & excep) {
      std::cerr << excep.what() << std::endl;
      std::cout << "{\"success\": false, \"error\": \"" << excep.what() << "\"}" << std::endl;
    }
  }

  dash::finalize();

  return EXIT_SUCCESS;
}
