#include <libdash.h>
#include <iostream>
#include <vector>

using std::cout;
using std::endl;
using std::cin;

using value       = struct{ uint row, col; };

using pattern_t   = dash::CSRPattern< 1, dash::ROW_MAJOR, long >;
using extent_t    = pattern_t::size_type;
using res_array_t = dash::Array<value, long, pattern_t>;

static int myid;


int main( int argc, char* argv[] )
{
  dash::init( &argc,&argv );
  myid = dash::myid( );

  
  std::vector<extent_t> local_sizes_source{20,30,20,30};

  pattern_t pattern_source( local_sizes_source );
  dash::Array<value, long, pattern_t> source ( pattern_source );

  int i = 0;
  
  for( value * src = source.lbegin(); src < source.lend(); ++src, ++i )
  {
    src->col = i * myid;
    src->row = i + myid;
  }
  
  
  value * target = new value[ 100 ];

  if(0==myid) cout << "before dash::copy" << endl;
  
  dash::copy(source.begin(), source.end(), target);
  
  if(0==myid) cout << "after dash::copy" << endl;
  
  
  dash::finalize( );
}