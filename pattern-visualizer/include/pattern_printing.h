#include <dash/Exception.h>

#include <input_param.h>
#include <pattern_creation.h>
#include <output_pattern.h>

template<TypeName          type,
         dash::dim_t       ndim,
         dash::MemArrange  arr>
void
make_and_print_constructor_pattern_h3(const json con_params,
                                   const cli_params output_params) {
  auto pattern = make_constructor_pattern<type,ndim,arr>(con_params);
  print_example(pattern, output_params, con_params);
}

template<TypeName          type,
         dash::dim_t       ndim>
void
make_and_print_constructor_pattern_h2(const dash::MemArrange arr_r,
                                      const json con_params,
                                      const cli_params output_params) {
  if(arr_r == dash::ROW_MAJOR) {
    make_and_print_constructor_pattern_h3<type,ndim,dash::ROW_MAJOR>(con_params,output_params);
  } else if(arr_r == dash::COL_MAJOR) {
    make_and_print_constructor_pattern_h3<type,ndim,dash::COL_MAJOR>(con_params,output_params);
  }
}

template<TypeName          type>
void
make_and_print_constructor_pattern_h1(const dash::dim_t ndim_r,
                                      const dash::MemArrange arr_r,
                                      const json con_params,
                                      const cli_params output_params) {
  if(ndim_r == 1) {
    make_and_print_constructor_pattern_h2<type,1>(arr_r,con_params,output_params);
  } else if(ndim_r == 2) {
    make_and_print_constructor_pattern_h2<type,2>(arr_r,con_params,output_params);
  } else if(ndim_r == 3) {
    make_and_print_constructor_pattern_h2<type,3>(arr_r,con_params,output_params);
  } else if(ndim_r == 4) {
    make_and_print_constructor_pattern_h2<type,4>(arr_r,con_params,output_params);
  } else if(ndim_r == 5) {
    make_and_print_constructor_pattern_h2<type,5>(arr_r,con_params,output_params);
  }
}

template<>
void
make_and_print_constructor_pattern_h1<csr>(const dash::dim_t ndim_r,
                                           const dash::MemArrange arr_r,
                                           const json con_params,
                                           const cli_params output_params) {
  if(ndim_r == 1) {
    make_and_print_constructor_pattern_h2<csr,1>(arr_r,con_params,output_params);
  } else {
    // CSRPattern is only implemented for NumDims == 1
    throw dash::exception::NotImplemented("CSRPattern is only implemented for "
                                          "NumDimensions == 1");
  }
}


void
make_and_print_constructor_pattern(const std::string type_r,
                                   const dash::dim_t ndim_r,
                                   const dash::MemArrange arr_r,
                                   const json con_params,
                                   const cli_params output_params) {
  if(type_r == "block") {
    make_and_print_constructor_pattern_h1<block>(ndim_r,arr_r,con_params,output_params);
  } else if(type_r == "tile") {
    make_and_print_constructor_pattern_h1<tile>(ndim_r,arr_r,con_params,output_params);
  } else if(type_r == "shiftTile") {
    make_and_print_constructor_pattern_h1<shiftTile>(ndim_r,arr_r,con_params,output_params);
  } else if(type_r == "seqTile") {
    make_and_print_constructor_pattern_h1<seqTile>(ndim_r,arr_r,con_params,output_params);
  /*} else if(type_r == "dynamic") {
    make_and_print_constructor_pattern<dynamic>(ndim_r,arr_r,con_params,output_params);*/
  } else if(type_r == "csr") {
    make_and_print_constructor_pattern_h1<csr>(ndim_r,arr_r,con_params,output_params);
  }
}

