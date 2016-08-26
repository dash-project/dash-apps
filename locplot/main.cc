
#include "json.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>


using json = nlohmann::json;
using std::cout;
using std::cerr;
using std::endl;


static std::map<std::string, std::string> scope_fills = {
#if 0
  { "GLOBAL", "#efefef" },
  { "NODE",   "#cfcfcf" },
  { "MODULE", "#e9c6af" },
  { "NUMA",   "#d89d96" },
  { "GROUP",  "#cf6262" },
  { "CORE",   "#87cdde" }
#else
  { "GLOBAL", "#ffffff" },
  { "NODE",   "#eff3f4" },
  { "MODULE", "#d0e0eb" },
  { "NUMA",   "#88abc2" },
  { "GROUP",  "#de8787" },
  { "CORE",   "#87cdde" }
#endif
};

static std::map<std::string, std::string> scope_strokes = {
#if 0
  { "GLOBAL", "#545454" },
  { "NODE",   "#232323" },
  { "MODULE", "#a05a2c" },
  { "NUMA",   "#944940" },
  { "GROUP",  "#800000" },
  { "CORE",   "#164450" }
#else
  { "GLOBAL", "#343434" },
  { "NODE",   "#343434" },
  { "MODULE", "#8aa0ab" },
  { "NUMA",   "#386583" },
  { "GROUP",  "#fe9a9a" },
  { "CORE",   "#287d93" }
#endif
};

static std::string fontstyle =
  " style=\"font-family:Lucida Console\" fill=\"#000000\" font-size=\"10\" ";

static std::string fontstyle_bold =
  " style=\"font-family:Lucida Console;font-weight:bold;\" fill=\"#000000\" font-size=\"10\" ";

typedef struct rect_s {
  int x;
  int y;
  int w;
  int h;
} rect_t;

typedef struct svg_node_s {
  rect_t      rect;
  std::string svg;
} svg_node_t;


std::string render_svg(
  const std::string & domain_tag,
  const json &        elem,
  int                 x,
  int                 y,
  int                 w,
  int                 h,
  int                 level)
{
  std::ostringstream os;

  std::string ind(level * 4, ' ');

  int tpad = 13;
  int fpad = 10;

  auto scope     = elem.find("scope");
  auto host      = elem.find("host");
  auto dom_level = elem.find("level");

  std::string hostname = "?";
  if (host != elem.end()) {
    hostname = *host;
  }

  int d_level = -1;
  if (dom_level != elem.end()) {
    d_level = *dom_level;
  }

  if (scope != elem.end()) {
    std::string scope_fill   = scope_fills[*scope];
    std::string scope_stroke = scope_strokes[*scope];
    std::string scope_name   = *scope;
    std::string rect_attr    = "";

    if (scope_name == "CORE") { scope_name = "UNIT"; }

    int col_0 = 60;
    int col_1 = 100;

    if (scope_name == "UNIT") {
      w     = 170;
      h     = 67;
      col_0 = 40;
      col_1 = 80;
    }
    if (scope_name == "GROUP") {
      rect_attr = " ry=\"8\" ";
    }
    os << ind << "<rect x=\"" << x << "\" y=\"" << y << "\""
              << " height=\"" << h << "\" width=\"" << w << "\""
              << rect_attr
              << " style=\""
                  << "fill:"   << scope_fill   << ";"
                  << "stroke:" << scope_stroke << ";"
                  << "stroke-width:1\" >"
              << "</rect>" << std::endl;

    os << ind << "<text"
              << " x=\"" << x + tpad << "\""
              << " y=\"" << y + tpad + fpad << "\""
              << fontstyle_bold
              << ">"
              << scope_name
              << "</text>" << endl;
    os << ind << "<text"
              << " x=\"" << x + tpad + col_0 << "\""
              << " y=\"" << y + tpad + fpad  << "\""
              << fontstyle
              << ">"
              << "L:" << d_level
              << "</text>" << endl;
    os << ind << "<text"
              << " x=\"" << x + tpad + col_1 << "\""
              << " y=\"" << y + tpad + fpad  << "\""
              << fontstyle
              << ">"
              << "[" << domain_tag << "]"
              << "</text>" << endl;

    if (scope_name == "NODE" || scope_name == "MODULE") {
      os << ind << "<text"
                << " x=\"" << x + tpad            << "\""
                << " y=\"" << y + tpad * 2 + fpad << "\""
                << fontstyle
                << ">"
                << elem["hwinfo"]["system_mb"]
                << " MB"
                << "</text>" << endl;
      os << ind << "<text"
                << " x=\"" << x + tpad + col_0    << "\""
                << " y=\"" << y + tpad * 2 + fpad << "\""
                << fontstyle
                << ">"
                << "host:" << hostname
                << "</text>" << endl;
    }
    if (scope_name == "NUMA") {
      os << ind << "<text"
                << " x=\"" << x + tpad            << "\""
                << " y=\"" << y + tpad * 2 + fpad << "\""
                << fontstyle
                << ">"
                << elem["hwinfo"]["numa_mb"]
                << " MB"
                << "</text>" << endl;
      os << ind << "<text"
                << " x=\"" << x + tpad + col_0    << "\""
                << " y=\"" << y + tpad * 2 + fpad << "\""
                << fontstyle
                << ">"
                << "id:" << elem["hwinfo"]["numa_id"]
                << "</text>" << endl;
    }
    if (scope_name == "UNIT") {
      int ncores   = elem["unit_locality"]["hwinfo"]["num_cores"];
      int nthreads = elem["unit_locality"]["hwinfo"]["threads"]["max"];
      os << ind << "<text"
                << " x=\"" << x + tpad << "\""
                << " y=\"" << y + tpad * 2 + fpad + 7 << "\""
                << fontstyle
                << ">"
                << "id:"
                << elem["unit_id"]["local_id"]
                << "</text>" << endl;
      os << ind << "<text"
                << " x=\"" << x + tpad << "\""
                << " y=\"" << y + tpad * 3 + fpad + 7 << "\""
                << fontstyle
                << ">"
                << ncores   << " x "
                << nthreads << " threads"
                << "</text>" << endl;
    }
  }

  return os.str();
}

svg_node_t render_domain(
  const std::string & domain_tag,
  const json &        elem,
  int                 x,
  int                 y,
  int                 level = 0)
{
  int w = 0;
  int h = 0;

  int pad  = 11;

  std::string ind(level * 4, ' ');
  std::ostringstream os;

  if (elem.is_object()) {
    auto domains  = elem.find("domains");
    auto scope    = elem.find("scope");
    bool vertical = true;
    if (scope != elem.end() && *scope == "NUMA") {
      vertical = false;
    }
    if (domains != elem.end()) {
      auto ndomains = std::distance(domains->begin(), domains->end());
      int  d_idx    = 0;
      for (auto d = domains->begin(); d != domains->end(); ++d) {
        int d_x = x;
        int d_y = y;
        if (vertical) {
          d_x += pad;
          d_y += h + 50;
        } else {
          d_x += w + pad;
          d_y += 50;
        }
        auto svg_node = render_domain(d.key(), d.value(),
                                      d_x, d_y,
                                      level+1);
        if (vertical) {
          h += svg_node.rect.h + (d_idx < ndomains - 1 ? pad : 0);
          w  = std::max(w, svg_node.rect.w);
        } else {
          w += svg_node.rect.w + (d_idx < ndomains - 1 ? pad : 0);
          h  = std::max(h, svg_node.rect.h);
        }
        os << svg_node.svg;
        d_idx++;
      }
      w += 2 * pad;
      h += 50 + pad;
    }
    else if (scope != elem.end() && *scope == "CORE") {
      w = 170;
      h = 67;
    }
  }

  rect_t     rect;
  svg_node_t node;
  rect.w    = w;
  rect.h    = h;
  node.rect = rect;
  node.svg  = render_svg(domain_tag, elem, x, y, rect.w, rect.h, level) +
              os.str();

  return node;
}

int main(int argc, char * argv[])
{
  if (argc < 2) {
    cerr << "no filename specified" << endl;
    return EXIT_FAILURE;
  }

  std::string   loc_json_file(argv[1]);
  std::ifstream loc_json_ifs(loc_json_file);
  if (!loc_json_ifs.is_open()) {
    cerr << "could not read file " << loc_json_file << endl;
    return EXIT_FAILURE;
  }
  std::string line;
  std::string loc_json_str;
  while (getline(loc_json_ifs, line)) {
    loc_json_str += line;
  }
  loc_json_ifs.close();

  auto loc_hierarchy = json::parse(loc_json_str);

  std::ostringstream os;
  os << "<svg xmlns=\"http://www.w3.org/2000/svg\"";
  os << " xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";

  os << render_domain(".", loc_hierarchy, 0, 0, 0).svg;

  cout << os.str() << endl;

  return EXIT_SUCCESS;
}
