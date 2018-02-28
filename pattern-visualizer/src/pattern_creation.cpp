#include <pattern_creation.h>
#include <main.h>
#include <input_param.h>

#include <dash/Pattern.h>
#include <dash/algorithm/SUMMA.h>

/*dash::TilePattern<2, dash::ROW_MAJOR, index_t>
make_summa_pattern(
  const cli_params                  & params,
  const dash::SizeSpec<2, extent_t> & sizespec,
  const dash::TeamSpec<2, index_t>  & teamspec)
{
  auto pattern = dash::make_pattern<
                   dash::summa_pattern_partitioning_constraints,
                   dash::summa_pattern_mapping_constraints,
                   dash::summa_pattern_layout_constraints >(
                     sizespec,
                     teamspec);

  if (params.tile_x >= 0 || params.tile_y >= 0) {
    // change tile sizes of deduced pattern:
    typedef decltype(pattern) pattern_t;
    pattern_t custom_pattern(sizespec,
                             dash::DistributionSpec<2>(
                               params.tile_x > 0
                               ? dash::TILE(params.tile_x)
                               : dash::NONE,
                               params.tile_y > 0
                               ? dash::TILE(params.tile_y)
                               : dash::NONE),
                             teamspec);
    pattern = custom_pattern;
  }
  return pattern;
}*/
