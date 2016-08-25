
#include "json.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using json = nlohmann::json;
using std::cout;
using std::cerr;
using std::endl;


std::string render_core_domain(
  const json & elem,
  int          x,
  int          y,
  int          w,
  int          h,
  int          level = 0)
{
  static std::string fontstyle =
    " style=\"font-family:Lucida Console\" fill=\"#000000\" font-size=\"10\" ";

  std::ostringstream os;
  std::string ind(level * 4, ' ');

  int pad  = 10;
  int tpad = 15;

  os << ind << "<rect x=\"" << x << "\" y=\"" << y << "\""
            << " height=\"" << h << "\" width=\"" << w << "\""
            << " style=\"fill:#ff5555;stroke:#124512;stroke-width:1\" >"
            << "</rect>" << std::endl;

  auto unit_id  = elem.find("unit_id");
  auto unit_loc = elem.find("unit_locality");

  int  nthreads = -1;
  if (unit_loc != elem.end()) {
//  nthreads = (*unit_loc)["threads"]["max"];
  }
  if (unit_id == elem.end()) {
    return "???";
  }
  int unit_local_id = (*unit_id)["local_id"];
  os << ind << "<text"
            << " x=\"" << x + tpad << "\""
            << " y=\"" << y + tpad << "\""
            << fontstyle
            << ">"
            << unit_local_id << " ("
            << " t:" << nthreads << ")"
            << "</text>" << endl;

  return os.str();
}

std::string render_domain(
  const json & elem,
  int          x,
  int          y,
  int          parent_w,
  int          parent_h,
  int          level = 0)
{
  std::ostringstream os;

  int w = parent_w;
  int h = parent_h;

  int pad  = 10;
  int tpad = 15;

  std::string ind(level * 4, ' ');

  static std::vector<std::string> scope_fills = {
    "#efefef", "#cfcfcf", "#afafaf", "#8f8f8f"
  };
  std::string scope_fill = scope_fills[
                             std::min<int>(level,
                                           scope_fills.size()-1)
                           ];

  static std::string fontstyle =
    " style=\"font-family:Lucida Console\" fill=\"#000000\" font-size=\"10\" ";

  os << ind << "<rect x=\"" << x << "\" y=\"" << y << "\""
            << " height=\"" << h << "\" width=\"" << w << "\""
            << " style=\"fill:" << scope_fill << ";"
            << " stroke:#005544;stroke-width:1\" >"
            << "</rect>" << std::endl;

  if (elem.is_object()) {
    auto scope   = elem.find("scope");
    auto domains = elem.find("domains");
    if (scope != elem.end()) {
      if (*scope == "CORE") {
        return render_core_domain(elem, x, y, w, h, level+1);
      }
      os << ind << "<text"
                << " x=\"" << x + tpad << "\""
                << " y=\"" << y + tpad << "\""
                << " fill=\"#000000\" font-size=\"10\"" << " >"
                << *scope
                << "</text>" << endl;
    }
    if (domains != elem.end()) {
      auto ndomains = std::distance(domains->begin(), domains->end());
      int  d_w      = (w / ndomains) - ((pad / ndomains) * (ndomains-1));
      int  d_h      = 40;
      int  d_idx    = 0;
      for (auto d = domains->begin(); d != domains->end(); ++d) {
        int d_x = x + (d_idx * (d_w + pad));
        int d_y = y + d_h + pad;

        os << render_domain(d.value(), d_x, d_y, d_w, d_h, level+1);

        os << ind << "<text"
                  << " x=\"" << d_x + tpad << "\""
                  << " y=\"" << d_y + tpad + 20 << "\""
                  << fontstyle
                  << ">"
                  << "domain:" << d.key()
                  << "</text>" << endl;

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

  os << render_domain(loc_hierarchy, 0, 0, 1200, 40);

  cout << os.str() << endl;

  return EXIT_SUCCESS;
}