/**
Dummy C++ app as a startpoint for writing dash applications
*/

#include <libdash.h>

int main(int argc, char** argv){
  // Initialize DASH
  dash::init(&argc, &argv);

  // do something

  // quit dash
  dash::finalize();
}
