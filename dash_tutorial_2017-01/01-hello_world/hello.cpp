#include <iostream>
#include <libdash.h>
using namespace std;

int main( int argc, char* argv[] ) {
  pid_t pid; char buf[100];

  dash::init( &argc, &argv );
  gethostname(buf, 100); pid = getpid();

  cout<<"'Hello world' from unit "<<
    dash::myid()<<" of "<<dash::size()<<
    " on "<<buf<<" pid="<<pid<<endl;

  dash::finalize();
}
