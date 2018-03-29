#include <libdash.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <thread>

using std::cout;
using std::endl;
using std::cin;

using value       = struct{ uint row, col; };

using pattern_t   = dash::CSRPattern< 1, dash::ROW_MAJOR, long >;
using extent_t    = pattern_t::size_type;
using res_array_t = dash::Array<value, long, pattern_t>;

static int myid;


using std::this_thread::sleep_for;
inline void __sleep( uint const mult = 10, uint const baseDur = 0 )
{
   uint SLEEP_TIME__ = myid * mult + baseDur;
    sleep_for(std::chrono::milliseconds(SLEEP_TIME__)); 
}



int main( int argc, char* argv[] )
{
  dash::init( &argc,&argv );
  myid = dash::myid( );

  
  std::vector<extent_t> local_sizes_source;
  
  int i;  
  
  for( i=0; i < dash::Team::All( ).size( ); ++i ){
    local_sizes_source.push_back(4+i);
  }

  pattern_t pattern_source( local_sizes_source );
  dash::Array<value, long, pattern_t> source ( pattern_source );

  i=0;  
  
  for( value * src = source.lbegin( ); src < source.lend( ); ++src, ++i )
  {
    src->col = i * myid;
    src->row = i + myid;
  }
  
  
  if(0==myid) cout << "size of array:" << source.size( ) << endl;
  value * target = new value[ source.size( ) ];
    
  
  // get Pattern
  auto patt = source.pattern( );
  
  __sleep(100);
  
  // iterate with GlobIt over Array Range
  for( auto GlobIt = source.begin( ); GlobIt < source.end( ); ++GlobIt ) 
  {
    auto localIndex = patt.local( GlobIt.pos( ) );
    cout << "id:" << myid
         << " globalIX:"           << std::setw(3) << GlobIt.pos( )
         << " local index/offset:" << std::setw(2) << localIndex.index
         << " on unit:"                            << localIndex.unit
         << "\n";
  }
  
  cout << endl;
  __sleep( );
  dash::barrier( );
  
  if(0==myid) cout << "before dash::copy" << endl;
  if(0==myid) dash::copy( source.begin( ), source.end( ), target );
  if(0==myid) cout <<  "after dash::copy" << endl;
  
  
  dash::finalize( );
}