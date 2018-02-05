#include <pattern_creation.h>
#include <main.h>
#include <input_param.h>

#include <dash/Pattern.h>
#include <dash/algorithm/SUMMA.h>

dash::ShiftTilePattern<2, dash::ROW_MAJOR, index_t>
make_shift_tile_pattern(
  const cli_params                  & params,
  const dash::SizeSpec<2, extent_t> & sizespec,
  const dash::TeamSpec<2, index_t>  & teamspec)
{
  // Example: -n 1680 1680 -u 28 1 -t 60 60
  typedef dash::ShiftTilePattern<2> pattern_t;
  pattern_t pattern(sizespec,
                    dash::DistributionSpec<2>(
                      dash::TILE(params.tile_y),
                      dash::TILE(params.tile_x)),
                    teamspec);
  return pattern;
}

dash::SeqTilePattern<2, dash::ROW_MAJOR, index_t>
make_seq_tile_pattern(
  const cli_params                  & params,
  const dash::SizeSpec<2, extent_t> & sizespec,
  const dash::TeamSpec<2, index_t>  & teamspec)
{
  // Example: -n 30 30 -u 4 1 -t 10 10
  typedef dash::SeqTilePattern<2> pattern_t;
  pattern_t pattern(sizespec,
                    dash::DistributionSpec<2>(
                      dash::TILE(params.tile_y),
                      dash::TILE(params.tile_x)),
                    teamspec);
  return pattern;
}

dash::TilePattern<2, dash::ROW_MAJOR, index_t>
make_tile_pattern(
  const cli_params                  & params,
  const dash::SizeSpec<2, extent_t> & sizespec,
  const dash::TeamSpec<2, index_t>  & teamspec)
{
  // Example: -n 30 30 -u 4 1 -t 10 10
  typedef dash::TilePattern<2> pattern_t;
  pattern_t pattern(sizespec,
                    dash::DistributionSpec<2>(
                      dash::TILE(params.tile_y),
                      dash::TILE(params.tile_x)),
                    teamspec);
  return pattern;
}

dash::Pattern<2, dash::ROW_MAJOR, index_t>
make_block_pattern(
  const cli_params                  & params,
  const dash::SizeSpec<2, extent_t> & sizespec,
  const dash::TeamSpec<2, index_t>  & teamspec)
{
  // Example: -n 30 30 -u 4 1 -t 10 10
  typedef dash::Pattern<2> pattern_t;
  pattern_t pattern(sizespec,
                    dash::DistributionSpec<2>(
                      (params.tile_y > 0
                       ? dash::BLOCKCYCLIC(params.tile_y)
                       : dash::NONE),
                      (params.tile_x > 0
                       ? dash::BLOCKCYCLIC(params.tile_x)
                       : dash::NONE)),
                    teamspec);
  return pattern;
}

dash::TilePattern<2, dash::ROW_MAJOR, index_t>
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
}
