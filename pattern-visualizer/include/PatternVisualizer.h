#ifndef PATTERN_VISUALIZER__pattern_visualizer_h
#define PATTERN_VISUALIZER__pattern_visualizer_h

#include <dash/dart/if/dart.h>

#include <sstream>
#include <array>
#include <iomanip>
#include <string>

namespace dash {
namespace tools {

/**
 *
 * Take a generic pattern instance and dump it as an JSON string.
 * To reduce size for higher-dimensional patterns it can be limited which
 * dimensions are included in the output.
 */
template<typename PatternT>
class PatternVisualizer
{
private:
  typedef PatternVisualizer<PatternT>      self_t;
  typedef typename PatternT::index_type    index_t;

private:
  const PatternT & _pattern;

public:
  /**
   * Constructs the Pattern Visualizer with a pattern instance.
   *
   * The pattern instance is constant. For a different pattern a new
   * PatternVisualizer has to be constructed.
   */
  PatternVisualizer(const PatternT & pat)
  : _pattern(pat)
  { }

  PatternVisualizer() = delete;
  PatternVisualizer(const self_t & other) = delete;

  /**
   * Outputs the pattern as a svg over the given output stream.
   * This method should only be called by a single unit.
   */
  void draw_pattern(
      std::ostream & os,
      /// Output the blocks of the pattern
      bool output_blocks = true,
      /// Output the memory layout of the pattern
      bool output_memlayout = true,
      // Only output the memlayout of the first unit
      bool memlayout_local = true,
      /// Reduce size by excluding extra dimensions.
      bool reduced = false,
      /// For higher dimensional patterns, defines which slice gets displayed
      std::array<index_t, PatternT::ndim()> coords = {},
      /// Defines which dimensions are written
      std::vector<index_t> dims = {}) {
    os << "{";

    // default dims to all dimensions
    if(dims.size() == 0) {
      dims.resize(PatternT::ndim());
      std::iota(dims.begin(),dims.end(),0);
    }

    // reduce to the first two dimensions specified
    if(reduced && dims.size() > 2) {
      dims.resize(2);
    }

    os << "\"dims\": [";
    for(auto it=dims.cbegin(); it != dims.cend(); ++it) {
      if(it != dims.cbegin()) {
        os << ",";
      }
      os << *it;
    }
    os << "]";

    if(output_blocks) {
      // default blocksize
      os << ",\n";
      os << "\"blocksize\": [";
      for(auto it=dims.cbegin(); it != dims.cend(); ++it) {
        if(it != dims.cbegin()) {
          os << ",";
        }
        os << _pattern.blocksize(*it);
      }
      os << "]";

      os << ",\n";
      draw_blocks(os, coords, dims);
    }

    // Redo in seperate function for irregular patterns
    // perhaps read information from memlayout
    //draw_tiles(os, sz, coords, dimx, dimy);

    if(output_memlayout) {
      os << ",\n";
      draw_memlayout(os, memlayout_local, coords, dims);
    }

    os << "}" << std::endl;
  }

  /**
   * Draws the separate tiles of the pattern
   */
  /*void draw_tiles(std::ostream & os,
                  const sizes & sz,
                  std::array<index_t, PatternT::ndim()> coords,
                  int dimx, int dimy) {
    for (int i = 0; i < _pattern.extent(dimx); i++) {
      for (int j = 0; j < _pattern.extent(dimy); j++) {

        os << "<rect x=\"" << (i * sz.gridx) << "\" y=\"" << (j * sz.gridy) << "\" ";
        os << "height=\"" << sz.tileszy << "\" width=\"" << sz.tileszx << "\" ";

        coords[dimx] = i;
        coords[dimy] = j;

        auto unit  = _pattern.unit_at(coords);
        auto loffs = _pattern.at(coords);

        os << tilestyle(unit);
        os << " tooltip=\"enable\" > ";
        os << " <title>Elem: (" << j << "," << i <<"),";
        os << " Unit " << unit;
        os << " Local offs. " << loffs;
        os << "</title>";

          //os << "<!-- i=" << i << " j=" << j << "--> ";
        os << "</rect>" << std::endl;
      }
    }
  }*/

  /**
   * Draws the blocks of the pattern
   */
  void draw_blocks(std::ostream & os,
                   std::array<index_t, PatternT::ndim()> coords,
                   std::vector<index_t> dims) {
    auto blockspec = get_blockspec(_pattern.blockspec());

    auto block_begin_coords = coords;
    auto block_begin_idx = _pattern.block_at(block_begin_coords);
    auto block_coords = blockspec.coords(block_begin_idx);

    os << "\"blocks\": [";

    bool first_iteration = true;
    for(auto rit=std::prev(dims.crend()); rit != dims.crend(); ++rit) {
      auto cur_dim = *rit;

      if(first_iteration) {
        block_coords[cur_dim] = 0;
        first_iteration = false;
      } else {
        ++block_coords[cur_dim];
        if(!(block_coords[cur_dim] < blockspec.extent(cur_dim))) {
          os << "]";
          continue;
        }
        os << "," << std::endl;
      }

      while(rit != dims.crbegin()) {
        --rit;
        cur_dim = *rit;

        block_coords[cur_dim] = 0;
        os << "[";
      }

      while(block_coords[cur_dim] < blockspec.extent(cur_dim)) {
        if(block_coords[cur_dim] != 0) {
          os << ",";
        }
        auto block_idx = blockspec.at(block_coords);
        auto block     = _pattern.block(block_idx);

        for(auto it=dims.cbegin(); it != dims.cend(); ++it) {
          block_begin_coords[*it] = block.offset(*it);
        }
        auto unit      = _pattern.unit_at(block_begin_coords);

        os << "{\"u\":" << unit;
        // include blocksize for underfilled blocks
        if(std::accumulate(dims.cbegin(),dims.cend(), false,
              [&](const bool & ret, const index_t & dim){
                return ret || (_pattern.blocksize(dim) != block.extent(dim));
              })) {
          os << ",\"s\":[";
          for(auto it=dims.cbegin(); it != dims.cend(); ++it) {
            if(it != dims.cbegin()) {
              os << ",";
            }
            os << block.extent(*it);
          }
          os << "]";
        }
        os << "}";

        ++block_coords[cur_dim];
      }

      os << "]";
    }
  }

  /**
   * Draws the memory layout for the current unit (usually unit 0)
   */
  void draw_memlayout(std::ostream & os,
                      // Only output for the first unit
                      bool local,
                      std::array<index_t, PatternT::ndim()> coords_slice,
                      std::vector<index_t> dims) {
    /*int startx, starty;
    int endx, endy;*/

    os << "\"memlayout\": [";

    auto num_units = _pattern.teamspec().size();
    for(dash::team_unit_t unit(0); unit < num_units; unit++) {
      if(unit != 0) {
        os << ",\n";
      }
      auto local_extents = _pattern.local_extents();//unit);
      auto local_size = std::accumulate(local_extents.begin(),
                                        local_extents.end(),
                                        1,std::multiplies<size_t>());
      bool contiguous = true;
      bool first_pass = true;
      os << "[";
      for(auto offset = 0; offset < local_size; offset++ ) {
        auto coords = _pattern.coords(_pattern.global(offset));//,unit));

        // compare if found coord is in current slice
        bool cur_slice = true;
        for(auto i=0; i < coords.size(); i++) {
          cur_slice = cur_slice &&
                      (coords_slice[i] == coords[i] ||
                       dims.end() != std::find(dims.begin(),dims.end(),i));
        }
        if(cur_slice) {
          if(!first_pass) {
            os << ",";
          }
          os << "{";
          // include offset if elements are skipped
          if(!contiguous) {
            os << "\"o\":" << offset << ",";
            contiguous = true;
          }
          os << "\"p\":[";
          for(auto it=dims.cbegin(); it != dims.cend(); ++it) {
            if(it != dims.cbegin()) {
              os << ",";
            }
            os << coords[*it];
          }
          os << "]}";
          first_pass = false;
        } else {
          contiguous = false;
        }
      }
      os << "]";

      // stop after unit 0
      if(local) {
        break;
      }
    }
    os << "]";
  }

private:

  template<typename BlockSpec_t>
  typename std::enable_if<
    BlockSpec_t::ndim::value >= 2,
    BlockSpec_t
  >::type
  get_blockspec(const BlockSpec_t & blockspec) {
    return blockspec;
  }

  template<typename BlockSpec_t>
  typename std::enable_if<
    BlockSpec_t::ndim::value == 1,
    CartesianIndexSpace<1, ROW_MAJOR, index_t>
  >::type
  get_blockspec(const BlockSpec_t & blockspec) {
    return CartesianIndexSpace<1, ROW_MAJOR, index_t>(blockspec.extents());
  }

};

} // namespace tools
} // namespace dash

#endif // PATTERN_VISUALIZER__pattern_visualizer_h
