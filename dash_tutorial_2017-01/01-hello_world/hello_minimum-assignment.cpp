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

    /* Assignment: Enlarge the distributed array such that every unit has the same number 
    of entries regardless of the number of units. Use values that result in a search time 
    of few seconds for the min_element() calls below for the slower case. */

    dash::Array<int> arr( /* fill me */ );

    /* Assignment: replace the following loop with a call to a STL or DASH template algorithm fill
    that performs the same thing. */

    for ( auto it= arr.lbegin(); it != arr.lend(); it++ ) {

        *it= dash::myid();
    }

    
    arr.barrier();
    start= std::chrono::system_clock::now();
    arr.barrier();
       
    /* Assignment: call the STL std::min_element() with arr's global iterator 
    on unit 0 just to show that one can. */
        
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

    /* Assignment: call the dash::min_element() counterpart as a collective operation 
    on all units. */

    arr.barrier();
    end= std::chrono::system_clock::now();
    
    if ( 0 == dash::myid() ) {
        cout << "dash::min_element() unit " << dash::myid() << " / " << dash::size() << " on " << buf << " pid= " << pid << " needed " <<
        std::chrono::duration_cast<std::chrono::milliseconds> (end-start).count() << " ms" << endl;
    }

    /* Assignment: Run with various numbers of units and watch how the runtimes behave. */
    
    dash::finalize();
}
