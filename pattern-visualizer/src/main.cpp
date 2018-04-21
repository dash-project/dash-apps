/**
 * \example ex.06.pattern-block-visualizer/main.cpp
 * Example demonstrating the instantiation of 
 * different patterns and their visualization.
 */ 

#include <main.h>
#include <input_param.h>
#include <pattern_printing.h>
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

        make_and_print_constructor_pattern(type_r,ndim_r,arr_r,con_params,output_params);
      }
    } catch (std::exception & excep) {
      std::cerr << excep.what() << std::endl;
      std::cout << "{\"success\": false, \"error\": \"" << excep.what() << "\"}" << std::endl;
    }
  }

  dash::finalize();

  return EXIT_SUCCESS;
}
