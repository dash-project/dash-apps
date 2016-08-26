
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
  { "GLOBAL", "#efefef" },
  { "NODE",   "#cfcfcf" },
  { "MODULE", "#e9c6af" },
  { "NUMA",   "#d89d96" },
  { "GROUP",  "#cf6262" },
  { "CORE",   "#87cdde" }
};

static std::map<std::string, std::string> scope_strokes = {
  { "GLOBAL", "#545454" },
  { "NODE",   "#232323" },
  { "MODULE", "#a05a2c" },
  { "NUMA",   "#944940" },
  { "GROUP",  "#800000" },
  { "CORE",   "#164450" }
};

static std::string fontstyle =
  " style=\"font-family:Lucida Console\" fill=\"#000000\" font-size=\"10\" ";

static std::string fontstyle_bold =
  " style=\"font-family:Lucida Console;font-weight:bold;\" fill=\"#000000\" font-size=\"10\" ";

std::string render_unit_domain(
  const std::string & domain_tag,
  const json &        elem,
  int                 x,
  int                 y,
  int                 parent_w,
  int                 parent_h,
  int                 level = 0)
{
  std::ostringstream os;

  int w = parent_w;
  int h = parent_h;

  int pad  = 10;
  int tpad = 12;

  std::string ind(level * 4, ' ');

  if (elem.is_object()) {
    auto scope   = elem.find("scope");
    auto host    = elem.find("host");
    auto domains = elem.find("domains");
  }

  return os.str();
}

std::string render_domain(
  const std::string & domain_tag,
  const json &        elem,
  int                 x,
  int                 y,
  int                 parent_w,
  int                 parent_h,
  int                 level = 0)
{
  std::ostringstream os;

  int w = parent_w;
  int h = parent_h;

  int pad  = 10;
  int tpad = 15;

  std::string ind(level * 4, ' ');

  if (elem.is_object()) {
    auto scope   = elem.find("scope");
    auto host    = elem.find("host");
    auto domains = elem.find("domains");
    std::string hostname = "?";
    if (host != elem.end()) {
      hostname = *host;
    }
    if (scope != elem.end()) {
      std::string scope_fill   = scope_fills[*scope];
      std::string scope_stroke = scope_strokes[*scope];
      std::string scope_name   = *scope;
      std::string rect_attr    = "";

      if (scope_name == "CORE") { scope_name = "UNIT"; }

      if (scope_name == "UNIT") {
        h         = 60;
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
                << " y=\"" << y + tpad << "\""
                << fontstyle_bold
                << ">"
                << scope_name
                << "</text>" << endl;
      os << ind << "<text"
                << " x=\"" << x + tpad + 40 << "\""
                << " y=\"" << y + tpad      << "\""
                << fontstyle
                << ">"
                << domain_tag
                << "</text>" << endl;

      if (scope_name == "NODE" || scope_name == "MODULE") {
        os << ind << "<text"
                  << " x=\"" << x + tpad + 90 << "\""
                  << " y=\"" << y + tpad << "\""
                  << fontstyle
                  << ">"
                  << "host:" << hostname
                  << "</text>" << endl;
        os << ind << "<text"
                  << " x=\"" << x + tpad     << "\""
                  << " y=\"" << y + tpad * 2 << "\""
                  << fontstyle
                  << ">"
                  << elem["hwinfo"]["system_mb"]
                  << " MB"
                  << "</text>" << endl;
      }
      if (scope_name == "NUMA") {
        os << ind << "<text"
                  << " x=\"" << x + tpad     << "\""
                  << " y=\"" << y + tpad * 2 << "\""
                  << fontstyle
                  << ">"
                  << elem["hwinfo"]["numa_mb"]
                  << " MB"
                  << "</text>" << endl;
      }
      if (scope_name == "UNIT") {
        int ncores   = elem["unit_locality"]["hwinfo"]["num_cores"];
        int nthreads = elem["unit_locality"]["hwinfo"]["threads"]["max"];
        os << ind << "<text"
                  << " x=\"" << x + tpad << "\""
                  << " y=\"" << y + tpad * 2 << "\""
                  << fontstyle
                  << ">"
                  << "id:"
                  << elem["unit_id"]["local_id"]
                  << "</text>" << endl;
        os << ind << "<text"
                  << " x=\"" << x + tpad << "\""
                  << " y=\"" << y + tpad * 3 << "\""
                  << fontstyle
                  << ">"
                  << ncores   << " x "
                  << nthreads << " threads"
                  << "</text>" << endl;
      }
    }
    if (domains != elem.end()) {
      auto ndomains = std::distance(domains->begin(), domains->end());
      int  d_w      = (w / ndomains) - ((pad / ndomains) * (ndomains-1));
      int  d_h      = 40;
      int  d_idx    = 0;
      for (auto d = domains->begin(); d != domains->end(); ++d) {
        int d_x = x + (d_idx * (d_w + pad));
        int d_y = y + d_h + pad;

        os << render_domain(d.key(), d.value(),
                            d_x, d_y, d_w, d_h, level+1);

        d_idx++;
      }
    }
  }

  return os.str();
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

  os << render_domain(".", loc_hierarchy, 0, 0, 4000, 40);

  cout << os.str() << endl;

  return EXIT_SUCCESS;
}
