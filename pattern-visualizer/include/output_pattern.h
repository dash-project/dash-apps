#ifndef PATTERN_VISUALIZER__output_pattern_h
#define PATTERN_VISUALIZER__output_pattern_h

#include <output_info.h>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include <PatternVisualizer.h>

/**
 * Create string describing of pattern instance.
 */
template<typename PatternType>
static
std::string pattern_to_string(
  const PatternType & pattern)
{
  typedef typename PatternType::index_type index_t;

  dash::dim_t ndim = pattern.ndim();

  std::string storage_order = pattern.memory_order() == dash::ROW_MAJOR
                              ? "ROW_MAJOR"
                              : "COL_MAJOR";

  std::array<index_t, 2> blocksize;
  blocksize[0] = pattern.blocksize(0);
  blocksize[1] = pattern.blocksize(1);

  std::ostringstream ss;
  ss << "dash::"
     << PatternType::PatternName
     << "<"
     << ndim << ","
     << storage_order << ","
     << typeid(index_t).name()
     << ">(" << '\n'
     << "        SizeSpec:  " << pattern.sizespec().extents()  << "," << '\n'
     << "        TeamSpec:  " << pattern.teamspec().extents()  << "," << '\n'
     << "        BlockSpec: " << pattern.blockspec().extents() << "," << '\n'
     << "        BlockSize: " << blocksize << " )";

  return ss.str();
}

/**
 * Create filename describing of pattern instance.
 */
template<typename PatternType>
static
std::string pattern_to_filename(
  const PatternType & pattern)
{
  typedef typename PatternType::index_type index_t;

  dash::dim_t ndim = pattern.ndim();

  std::string storage_order = pattern.memory_order() == dash::ROW_MAJOR
                              ? "ROW_MAJOR"
                              : "COL_MAJOR";

  auto sspc = pattern.sizespec();
  auto tspc = pattern.teamspec();
  auto bspc = pattern.blockspec();

  std::ostringstream ss;
  ss << PatternType::PatternName
     << "--"
     << ndim << "-"
     << storage_order << "-"
     << typeid(index_t).name()
     << "--"
     << "size-"   << sspc.extent(0) << "x" << sspc.extent(1) << "--"
     << "team-"   << tspc.extent(0) << "x" << tspc.extent(1) << "--"
     << "blocks-" << bspc.extent(0) << "x" << bspc.extent(1)
     << ".svg";

  return ss.str();
}

template<typename PatternT>
static
void print_example(
  const PatternT   & pattern,
  const cli_params & output_params,
  const json params)
{
  typedef typename PatternT::index_type index_t;

  auto pattern_file = pattern_to_filename(pattern);
  auto pattern_desc = pattern_to_string(pattern);
  //print_pattern_metrics(pattern);

  dash::tools::PatternVisualizer<PatternT> pv(pattern);
  // pv.set_title(pattern_desc);

  std::cerr << "Generating visualization of "
            << '\n'
            << "    " << pattern_desc
            << std::endl;

  std::cout << "{\"success\": true, \"name\": \"" << pattern_file << "\", ";
  std::cout <<  "\"pattern\": ";
  pv.draw_pattern(std::cout);
  std::cout << "}" << std::endl;
}

#endif // PATTERN_VISUALIZER__output_pattern_h
