/**
Dummy C++ app as a startpoint for writing dash applications
*/
#include <iostream>

#include <libdash.h>

int main(int argc, char** argv){
  // Initialize DASH
  dash::init(&argc, &argv);

  // do something
  std::cout << "I'm Unit " << dash::myid() << std::endl;  

  // quit dash
  dash::finalize();
}
