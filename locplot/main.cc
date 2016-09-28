
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
  { "GLOBAL",  "#ffffff" },
  { "NODE",    "#eff3f4" },
  { "MODULE",  "#d0e0eb" },
  { "PACKAGE", "#d0e0eb" },
  { "NUMA",    "#c6e9af" },
  { "CACHE",   "#c7dae0" },
  { "GROUP",   "#de8787" },
  { "CORE",    "#87cdde" }
};

static std::map<std::string, std::string> scope_strokes = {
  { "GLOBAL",  "#343434" },
  { "NODE",    "#343434" },
  { "MODULE",  "#8aa0ab" },
  { "PACKAGE", "#8aa0ab" },
  { "NUMA",    "#5aa02c" },
  { "CACHE",   "#83b4c0" },
  { "GROUP",   "#fe9a9a" },
  { "CORE",    "#287d93" }
};

static std::string fontstyle_regular =
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

std::string svg_text(
  const std::string & text,
  int                 x,
  int                 y,
  const std::string & style = fontstyle_regular)
{
  std::ostringstream os;
  os << "<text"
     << " x=\"" << x << "\""
     << " y=\"" << y << "\""
     << style
     << ">"
     << text
     << "</text>" << endl;
  return os.str();
}

template <typename T>
std::string svg_text(
  const std::string & label,
  T                   value,
  int                 x,
  int                 y,
  const std::string & style = fontstyle_regular)
{
  std::ostringstream os;
  os << "<text"
     << " x=\"" << x << "\""
     << " y=\"" << y << "\""
     << style
     << ">"
     << label << ":" << value
     << "</text>" << endl;
  return os.str();
}

std::string svg_rect(
  int x,
  int y,
  int w,
  int h,
  const std::string & fill_color,
  const std::string & stroke_color,
  const std::string & attr = "")
{
  std::ostringstream os;
  os << "<rect x=\"" << x << "\" y=\"" << y << "\""
     << " height=\"" << h << "\" width=\"" << w << "\""
     << attr
     << " style=\""
         << "fill:"   << fill_color   << ";"
         << "stroke:" << stroke_color << ";"
         << "stroke-width:1\" >"
     << "</rect>" << std::endl;
  return os.str();
}

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
    int col_1 = 130;

    if (scope_name == "UNIT") {
      w     = 170;
      h     = 77;
    }
    if (scope_name == "GROUP") {
      rect_attr = " ry=\"8\" ";
    }
    os << ind << svg_rect(x, y, w, h, scope_fill, scope_stroke, rect_attr);
    os << ind << svg_text(scope_name,
                          x + tpad, y + tpad + fpad, fontstyle_bold);
    os << ind << svg_text("L", d_level,
                          x + tpad + col_0, y + tpad + fpad);
    os << ind << svg_text(std::string("[") + domain_tag + "]",
                          x + tpad, y + (tpad * 2) + fpad);

    auto elem_hwinfo = elem.find("hwinfo");
    std::ostringstream shared_mem_kb;
    if (elem_hwinfo != elem.end()) {
      auto hwinfo_shared_mem_kb = elem_hwinfo->find("shmem");
      if (hwinfo_shared_mem_kb != elem_hwinfo->end()) {
        shared_mem_kb << *hwinfo_shared_mem_kb << " KB";
      }
    }

    if (scope_name == "NODE" || scope_name == "MODULE") {
      std::ostringstream system_mb;
      system_mb << shared_mem_kb.str();
      os << ind << svg_text(std::string("host:") + hostname,
                            x + tpad + col_1, y + tpad + fpad);
      os << ind << svg_text(system_mb.str(),
                            x + tpad, y + tpad * 3 + fpad);
    }
    if (scope_name == "NUMA") {
      std::ostringstream numa_mb;
      numa_mb << shared_mem_kb.str();
      std::string numa_id = "?";
      if (elem_hwinfo != elem.end() &&
          elem_hwinfo->find("numa_id") != elem_hwinfo->end()) {
        numa_id = (*elem_hwinfo)["numa_id"];
      }
      os << ind << svg_text("id", numa_id,
                            x + tpad + col_1, y + tpad + fpad);
      os << ind << svg_text(numa_mb.str(),
                            x + tpad, y + tpad * 3 + fpad);
    }
    if (scope_name == "CACHE") {
      std::ostringstream cache_size;
      cache_size << shared_mem_kb.str();
      os << ind << svg_text(cache_size.str(),
                            x + tpad, y + tpad * 3 + fpad);
    }
    if (scope_name == "UNIT") {
      int ncores   = elem["unit_loc"]["hwinfo"]["num_cores"];
      int nthreads = elem["unit_loc"]["hwinfo"]["threads"]["max"];
      int cpu_id   = elem["unit_loc"]["hwinfo"]["cpu_id"];
      int unit_id  = elem["unit_id"]["local_id"];

      os << ind << svg_text("id", unit_id,
                            x + tpad, y + tpad * 3 + fpad + 7);
      os << ind << svg_text("CPU", cpu_id,
                            x + tpad + 60, y + tpad * 3 + fpad + 7);
      os << ind << svg_text("cores", ncores,
                            x + tpad, y + tpad * 4 + fpad + 7);
      os << ind << svg_text("t/c", nthreads,
                            x + tpad + 60, y + tpad * 4 + fpad + 7);
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
    bool nested   = true;
    if (domains != elem.end()) {
      auto ndomains = std::distance(domains->begin(), domains->end());
      if (scope != elem.end() &&
          (*scope == "NUMA" || *scope == "CACHE" || *scope == "PACKAGE")) {
        vertical = false;
      }
      if (scope != elem.end() && (*scope == "CACHE") && ndomains > 0) {
        nested = false;
      }
      int  d_idx = 0;
      for (auto d = domains->begin(); d != domains->end(); ++d) {
        int d_x = x;
        int d_y = y;

        if (nested) {
          if (vertical) {
            d_x += pad;
            d_y += h + 60;
          } else {
            d_x += w + pad;
            d_y += 60;
          }
        } else {
          // draw stacked domains, horizontal
          d_x += w;
          d_y += 60;
        }

        auto svg_node = render_domain(d.key(), d.value(),
                                      d_x, d_y,
                                      level+1);
        if (nested) {
          // draw nested domains
          if (vertical) {
            h += svg_node.rect.h + (d_idx < ndomains - 1 ? pad : 0);
            w  = std::max(w, svg_node.rect.w);
          } else {
            w += svg_node.rect.w + (d_idx < ndomains - 1 ? pad : 0);
            h  = std::max(h, svg_node.rect.h);
          }
        } else {
          // draw stacked domains, horizontal
          w += svg_node.rect.w + (d_idx < ndomains - 1 ? pad : 0);
          h  = std::max(h, svg_node.rect.h);
        }
        os << svg_node.svg;
        d_idx++;
      }
      if (nested) {
        w += 2 * pad;
        h += 60 + pad;
      } else {
        h += 60;
      }
    }
    else if (scope != elem.end() && *scope == "CORE") {
      w = 170;
      h = 77;
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
