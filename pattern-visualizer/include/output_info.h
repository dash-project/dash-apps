#ifndef PATTERN_VISUALIZER__output_info_h
#define PATTERN_VISUALIZER__output_info_h

#include <sstream>

#include <dash/util/PatternMetrics.h>
#include <dash/pattern/internal/PatternLogging.h>

template<typename PatternT>
static
void print_pattern_metrics(const PatternT & pattern)
{
  dash::util::PatternMetrics<PatternT> pm(pattern);

  size_t block_kbytes = pattern.blocksize(0) * pattern.blocksize(1) *
                        sizeof(double) /
                        1024;

  std::cerr << "Pattern Metrics:"
            << '\n'
            << "    Partitioning:"
            << '\n'
            << "        block size:         " << block_kbytes << " KB"
            << '\n'
            << "        number of blocks:   " << pm.num_blocks()
            << '\n'
            << "    Mapping imbalance:"
            << '\n'
            << "        min. blocks/unit:   " << pm.min_blocks_per_unit()
            << " = " << pm.min_elements_per_unit() << " elements"
            << '\n'
            << "        max. blocks/unit:   " << pm.max_blocks_per_unit()
            << " = " << pm.max_elements_per_unit() << " elements"
            << '\n'
            << "        imbalance factor:   " << std::setprecision(4)
                                              << pm.imbalance_factor()
            << '\n'
            << "        balanced units:     " << pm.num_balanced_units()
            << '\n'
            << "        imbalanced units:   " << pm.num_imbalanced_units()
            << '\n'
            << std::endl;
}

template<typename PatternT>
static
void log_pattern_mapping(const PatternT & pattern)
{
  typedef PatternT pattern_t;

  dash::internal::print_pattern_mapping(
    "pattern.unit_at", pattern, 3,
    [](const pattern_t & _pattern, int _x, int _y) -> dart_unit_t {
      return _pattern.unit_at(std::array<index_t, 2> {{_x, _y}});
    });
  dash::internal::print_pattern_mapping(
    "pattern.global_at", pattern, 3,
    [](const pattern_t & _pattern, int _x, int _y) -> index_t {
      return _pattern.global_at(std::array<index_t, 2> {{_x, _y}});
    });
  dash::internal::print_pattern_mapping(
    "pattern.local", pattern, 10,
    [](const pattern_t & _pattern, int _x, int _y) -> std::string {
      auto lpos = _pattern.local(std::array<index_t, 2> {{_x, _y}});
        std::ostringstream ss;
        ss << lpos.unit << ":" << lpos.coords;
        return ss.str();
    });
  dash::internal::print_pattern_mapping(
    "pattern.at", pattern, 3,
    [](const pattern_t & _pattern, int _x, int _y) -> index_t {
      return _pattern.at(std::array<index_t, 2> {{_x, _y}});
    });
  dash::internal::print_pattern_mapping(
    "pattern.block_at", pattern, 3,
    [](const pattern_t & _pattern, int _x, int _y) -> index_t {
      return _pattern.block_at(std::array<index_t, 2> {{_x, _y}});
    });
  dash::internal::print_pattern_mapping(
    "pattern.block.offset", pattern, 5,
    [](const pattern_t & _pattern, int _x, int _y) -> std::string {
      auto block_idx = _pattern.block_at(std::array<index_t, 2> {{_x, _y}});
        auto block_vs  = _pattern.block(block_idx);
        std::ostringstream ss;
        ss << block_vs.offset(0) << "," << block_vs.offset(1);
        return ss.str();
    });
  dash::internal::print_pattern_mapping(
    "pattern.local_index", pattern, 3,
    [](const pattern_t & _pattern, int _x, int _y) -> index_t {
      return _pattern.local_index(std::array<index_t, 2> {{_x, _y}}).index;
    });
}

#endif // PATTERN_VISUALIZER__output_info_h
