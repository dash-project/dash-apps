#ifndef PATTERN_VISUALIZER__pattern_creation_h
#define PATTERN_VISUALIZER__pattern_creation_h

#include <main.h>
#include <input_param.h>

#include <dash/Pattern.h>

dash::ShiftTilePattern<2, dash::ROW_MAJOR, index_t>
make_shift_tile_pattern(
  const cli_params                  & params,
  const dash::SizeSpec<2, extent_t> & sizespec,
  const dash::TeamSpec<2, index_t>  & teamspec);

dash::SeqTilePattern<2, dash::ROW_MAJOR, index_t>
make_seq_tile_pattern(
  const cli_params                  & params,
  const dash::SizeSpec<2, extent_t> & sizespec,
  const dash::TeamSpec<2, index_t>  & teamspec);

dash::TilePattern<2, dash::ROW_MAJOR, index_t>
make_tile_pattern(
  const cli_params                  & params,
  const dash::SizeSpec<2, extent_t> & sizespec,
  const dash::TeamSpec<2, index_t>  & teamspec);

dash::Pattern<2, dash::ROW_MAJOR, index_t>
make_block_pattern(
  const cli_params                  & params,
  const dash::SizeSpec<2, extent_t> & sizespec,
  const dash::TeamSpec<2, index_t>  & teamspec);

dash::TilePattern<2, dash::ROW_MAJOR, index_t>
make_summa_pattern(
  const cli_params                  & params,
  const dash::SizeSpec<2, extent_t> & sizespec,
  const dash::TeamSpec<2, index_t>  & teamspec);

#endif // PATTERN_VISUALIZER__pattern_creation_h
