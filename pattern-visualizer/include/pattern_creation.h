#ifndef PATTERN_VISUALIZER__pattern_creation_h
#define PATTERN_VISUALIZER__pattern_creation_h

#include <main.h>
#include <input_param.h>

#include <dash/Cartesian.h>
#include <dash/Dimensional.h>
#include <dash/TeamSpec.h>
#include <dash/Pattern.h>
//#include <dash/pattern/DynamicPattern.h>

enum TypeName { block, tile, shiftTile, seqTile, /*dynamic,*/ csr };

template<TypeName          type,
         dash::dim_t       ndim,
         dash::MemArrange  arr>
struct PatternType {};

template<dash::dim_t      ndim,
         dash::MemArrange arr>
struct PatternType<block,ndim,arr> {
  using type = dash::BlockPattern<ndim,arr>;
};

template<dash::dim_t      ndim,
         dash::MemArrange arr>
struct PatternType<tile,ndim,arr> {
  using type = dash::TilePattern<ndim,arr>;
};

template<dash::dim_t      ndim,
         dash::MemArrange arr>
struct PatternType<shiftTile,ndim,arr> {
  using type = dash::ShiftTilePattern<ndim,arr>;
};

template<dash::dim_t      ndim,
         dash::MemArrange arr>
struct PatternType<seqTile,ndim,arr> {
  using type = dash::SeqTilePattern<ndim,arr>;
};

/*template<dash::dim_t      ndim,
         dash::MemArrange arr>
struct PatternType<dynamic,ndim,arr> {
  using type = dash::DynamicPattern<ndim,arr>;
};*/

template<dash::dim_t      ndim,
         dash::MemArrange arr>
struct PatternType<csr,ndim,arr> {
  using type = dash::CSRPattern<ndim,arr>;
};

template<dash::dim_t ndim>
dash::SizeSpec<ndim>
make_size_spec(const json params) {
  std::array<extent_t, ndim> extents = {};
  std::transform(params.begin(),params.end(),extents.begin(),
                 [](const json & value) {
                   const std::string val = value;
                   return std::stoi(val);
                 });
  return dash::SizeSpec<ndim>(extents);
}

template<dash::dim_t ndim>
dash::DistributionSpec<ndim>
make_dist_spec(const json params) {
  std::array<dash::Distribution, ndim> dists;
  std::transform(params.begin(),params.end(),dists.begin(),
                 [](const json & value) {
                   std::string size;
                   const std::string name = value.at("name");
                   if(name == "BLOCKED") {
                     return dash::BLOCKED;
                   } else if(name == "CYCLIC") {
                     return dash::CYCLIC;
                   } else if(name == "NONE") {
                     return dash::NONE;
                   } else if(name == "TILE") {
                     size = value.at("params")[0].at("value");
                     return dash::TILE(std::stoi(size));
                   } else if(name == "BLOCKCYCLIC") {
                     size = value.at("params")[0].at("value");
                     return dash::BLOCKCYCLIC(std::stoi(size));
                   }
                 });
  return dash::DistributionSpec<ndim>(dists);
}

template<dash::dim_t ndim>
dash::TeamSpec<ndim>
make_team_spec(const json params) {
  std::array<typename dash::TeamSpec<ndim>::size_type,ndim> extents;
  std::fill(extents.begin(),extents.end(),1);
  const std::string name = params.at("value").at("name");
  if(name == "yes") {
    const std::string num_units = params.at("value").at("params")[0].at("value");
    extents[0] = std::stoi(num_units);
    dash::TeamSpec<ndim> team_spec(extents);
    team_spec.balance_extents();
    return team_spec;
  } else if(name == "no") {
    auto team_params = params.at("value").at("params")[0].at("value");
    std::transform(team_params.begin(),team_params.end(),extents.begin(),
                   [](const json & value) {
                     const std::string val = value;
                     return std::stoi(val);
                   });
    dash::TeamSpec<ndim> team_spec(extents);
    return team_spec;
  }
}


template<TypeName          type,
         dash::dim_t       ndim,
         dash::MemArrange  arr>
typename PatternType<type,ndim,arr>::type
make_constructor_pattern(json params) {
  auto extent_params = getParamFromGroup(params,"ext_group").at("content");
  auto size_spec = make_size_spec<ndim>(getParamFromGroup(extent_params,"size").at("value"));
  auto dist_spec = make_dist_spec<ndim>(getParamFromGroup(extent_params,"dist").at("value"));
  auto team_spec = make_team_spec<ndim>(getParamFromGroup(params,"balance_team"));

  typename PatternType<type,ndim,arr>::type pattern(size_spec,dist_spec,team_spec);

  return pattern;
}


//make_algorithm_pattern();

/*dash::TilePattern<2, dash::ROW_MAJOR, index_t>
make_summa_pattern(
  const cli_params                  & params,
  const dash::SizeSpec<2, extent_t> & sizespec,
  const dash::TeamSpec<2, index_t>  & teamspec);*/

#endif // PATTERN_VISUALIZER__pattern_creation_h
