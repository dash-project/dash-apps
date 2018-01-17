#include <iostream>
#include <libdash.h>
using namespace std;

int main( int argc, char* argv[] ) {

  dash::init( &argc, &argv );

  dash::Array<int> arr(100);

  if ( 1 == dash::myid() ) {
    for ( auto i= 0; i < arr.size(); i++ )
      arr[i]= i;
  }

  dash::barrier();

  if ( 0 == dash::myid() ) {
    for( auto el: arr ) {
      cout << (int)el << " ";
    }
    cout << endl;
  }

  dash::barrier();

  for ( auto it= arr.lbegin(); it != arr.lend(); it++ ) {
    *it= dash::myid();
  }

  arr.barrier();

  if ( 0 == dash::myid() ) {
    for ( auto el: arr ) {
      cout<<(int)el<<" ";
    }
    cout<<endl;
  }

  dash::finalize();
}
