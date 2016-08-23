
#include "json.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using json = nlohmann::json;
using std::cout;
using std::cerr;
using std::endl;


std::string render_domain(
  const json & elem,
  int          parent_w,
  int          parent_h)
{
  std::ostringstream os;

  int w = parent_w;
  int h = parent_h;

  int x = 5;
  int y = 5;

  os << "<rect x=\"" << x << "\" y=\"" << y << "\" "
     << "height=\"" << h << "\" width=\"" << w << "\" >"
     << endl;

  for (auto it = elem.begin(); it != elem.end(); ++it) {
    if (it->is_string()) {
      os << *it << endl;
    }
    if (it->is_object()) {
      if (it.key() == "NODE") {
        os << "NODE: " << it.value() << endl;
      } else {
        os << render_domain(*it, w, h);
      }
    }
  }

  os << endl;
  os << "</rect>" << std::endl;

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

  for (auto & elem : loc_hierarchy) {
    os << render_domain(elem, 100, 100);
  }

  cout << os.str() << endl;

  return EXIT_SUCCESS;
}
