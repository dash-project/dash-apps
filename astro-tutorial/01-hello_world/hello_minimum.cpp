#include <iostream>
#include <libdash.h>
#include <chrono>

using namespace std;

int main( int argc, char* argv[] ) {

    pid_t pid;
    char buf[100];
    std::chrono::time_point<std::chrono::system_clock> start, end;

    dash::init( &argc, &argv );
    gethostname( buf, 100 );
    pid = getpid();

    dash::Array<int> arr( 1000000 * dash::size() );

    /*
    for ( auto it= arr.lbegin(); it != arr.lend(); it++ ) {

        *it= dash::myid();
    }
    */

    // std::fill( arr.local.begin(), arr.local.end(), (int) dash::myid() );
    // std::fill( arr.lbegin(), arr.lend(), (int) dash::myid() );
    dash::fill( arr.begin(), arr.end(), (int) dash::myid() );

    arr.barrier();
    start= std::chrono::system_clock::now();
    arr.barrier();
    if ( 0 == dash::myid() ) {

        std::min_element( arr.begin(), arr.end() );
    }
    arr.barrier();
    end= std::chrono::system_clock::now();

    if ( 0 == dash::myid() ) {
        cout << "std::min_element() unit " << dash::myid() << " / " << dash::size() << 
        " on " << buf << " pid= " << pid << " needed " <<
        std::chrono::duration_cast<std::chrono::milliseconds> (end-start).count() << " ms" << endl;
    }

    arr.barrier();
    start= std::chrono::system_clock::now();
    arr.barrier();
    dash::min_element( arr.begin(), arr.end() );
    arr.barrier();
    end= std::chrono::system_clock::now();

    if ( 0 == dash::myid() ) {
        cout << "dash::min_element() unit " << dash::myid() << " / " << dash::size() << " on " << buf << " pid= " << pid << " needed " <<
        std::chrono::duration_cast<std::chrono::milliseconds> (end-start).count() << " ms" << endl;
    }

    dash::finalize();
}
