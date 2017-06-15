#ifndef DASH_ENABLE_LOGGING
#define DASH_ENABLE_LOGGING
#endif

#include <libdash.h>
#include <chrono>
#include <thread>
#include "../Terminal_Color.h"

//#define DEBUG
#define SLEEP_TIME            30 //that's the sleep time before DEBUG IO
#define MAX_KEY               99
#define MIN_NUM_ELEM_PER_UNIT 10


using std::cout;
using std::endl;
using std::vector;

#ifdef DEBUG
  using std::this_thread::sleep_for;
#endif

using uint       = unsigned  int;
using uchar      = unsigned char;

using MTRX_TYPE  = uchar;

using Point      = struct{ MTRX_TYPE value; uint row, col;};
using pointRange = struct{ Point * curr, * end;           };
using unitRange  = struct{ MTRX_TYPE begin, end;          };


typedef dash::CSRPattern< 1, dash::ROW_MAJOR, int > pattern_t;
typedef typename pattern_t::size_type extent_t;

template< typename T >
void print2d( T& mat ) {
  for( int i = 0; i < mat.extent(0); i++ ) {
    for( int j = 0; j < mat.extent(1); j++ ) {
      cout << std::setw(3) << static_cast<uint>( mat(i,j) )<< " ";
    }
    cout << endl;
  }
}

bool operator<(const Point& lhs, const Point& rhs)
{
  return lhs.value < rhs.value;
}

std::ostream& operator<<(std::ostream& os, const Point& p)
{
  #ifdef DEBUG
    return os << "(" << fmt( p.value, FCYN ) << "," << fmt( p.row, FGREEN ) << "," << fmt( p.col, FGREEN ) << ")-";
  #else
    return os << p.row << " " << p.col << endl;
  #endif
}


//template<typename T = MTRX_TYPE>
inline void winnow(uint nrows, uint ncols, uint nelts, int myid){

  dash::Team & team = dash::Team::All();
  size_t nUnits     = team.size();
  
  // create global histo array and 
  dash::Array<uint> histo( (MAX_KEY + 1) * nUnits, dash::BLOCKED );

  vector<Point > pointsLocal;
  vector<size_t> foundLclCpy( nUnits );
  
  
  
  /* read the matrix and boolean matrix from std in and generate a histogram and 
   * how much elements each unit found
   */
  size_t foundAllSize = readMatricesFromStdIn( nrows, ncols, pointsLocal, foundLclCpy, histo, nUnits, myid );

  
  // in this array will be the distribution info for creating buckets
  vector<unitRange> distr( nUnits );
  
  /* calculate global histogram and thereof a distribution pattern */
  calcGlobDist( histo, distr, foundAllSize, nUnits, myid, team );
  
  /* each unit creates buckets that are send to the responsible unit later */
  //createBuckets( myid,team )
  
  
}


template< typename T = MTRX_TYPE >
inline size_t readMatricesFromStdIn (
                 uint const   nrows,
                 uint const   ncols,
        vector<Point>       & pointsLocal,
       vector<size_t>       & foundLclCpy, 
  dash::Array< uint >       & histo,
               size_t         nUnits,
                  int const   myid
){
    dash::NArray<T   , 2> matrix( nrows, ncols );
    dash::NArray<bool, 2> mask  ( nrows, ncols );
    dash::Array <size_t > found ( nUnits );    
    
    // read input matrix from stdin
    if( 0 == myid ){
      T tmp;
      for ( auto i : matrix ){
        scanf( "%d", &tmp )  , i = tmp;
      }
      bool tmpB;
      for ( auto i : mask ){
        scanf( "%d", &tmpB ) , i = tmpB;
      }
    }
    
    dash::barrier( );
    
    // initialize histo with 0
    for( uint * i = histo.lbegin(); i < histo.lend(); ++i) {
      *i = 0;
    }
    

    #ifdef DEBUG  // print error message if mask's and matrix's local size aren't identical
      if( mask.lend( ) - mask.lbegin( ) != matrix.lend( ) - matrix.lbegin( )) {
        cout << "On unit " << myid << " the local sizes of matrix and mask differ!\naborted on this unit" << endl;
        return -1;
      }
    #endif
    
    
    auto globIndex = matrix.pattern( ).global( {0,0} );
    uint gRow      = globIndex[0];
    uint gCol      = globIndex[1];
    T * matrEl     = matrix.lbegin( );
    
    // read in local part of mask - matrix combination
    for ( bool * maskEl = mask.lbegin( );  maskEl < mask.lend( );  ++maskEl , ++matrEl, ++gCol ) {
      if( gCol == ncols ) gCol = 0, ++gRow;
      
      if( *maskEl       ) {
        pointsLocal.push_back( Point{ *matrEl, gRow, gCol } );
        ++histo.local[*matrEl];
      }
    }
    
    found[myid] = pointsLocal.size();  

  
    // wait for all to set found
  dash::barrier( );
  dash::copy(found.begin( ), found.end( ), foundLclCpy.data( ) );
  
  
  size_t foundAllSize = 0;
  for( size_t * i = foundLclCpy.data( ); i < foundLclCpy.data( ) + nUnits; ++i ) {
    foundAllSize += *i;
  }

  #ifdef DEBUG  // print points found local after sort
    sleep_for(std::chrono::milliseconds(SLEEP_TIME)); 
    
    cout << "#" << fmt( myid, FBLUE, 2 ) << ": ";
    for( auto it : pointsLocal) {
      cout << it;
    } cout << endl;

    sleep_for(std::chrono::milliseconds(SLEEP_TIME)); 
    cout << "#" << fmt( myid, FBLUE, 2 ) << ": found:" << fmt( pointsLocal.size(),FRED , 2 )  << endl;
  
    sleep_for(std::chrono::milliseconds(SLEEP_TIME)); 
    if( 0 == myid ) cout << "foundAllSize: " << fmt( foundAllSize, FRED ) << endl;
  #endif
  
  return foundAllSize;
}


inline void calcGlobDist(
  dash::Array< uint > & histo,
    vector<unitRange> & distr,
               size_t   foundAllSize,
               size_t   nUnits,
                  int   myid,
           dash::Team & team
){
  
  if( 0 != myid ) {
    dash::transform<uint>(
      histo.lbegin     ( ) ,
      histo.lend       ( ) ,
      histo.begin      ( ) ,
      histo.begin      ( ) ,
      dash::plus<uint> ( )
    );
  }
  
  
  // unit 0 have to wait for rest to finish adding their values
  dash::barrier();
  
  if( 0 == myid ) { // calculate bucket distribution

    /* calculate how much elements each unit should hold ideally
     * increment by one for safety garuantees in distribution
     * that is the assumption (ideal * nUnits > foundAllSize) == true
     */
    uint ideal = std::max( static_cast<size_t>( MIN_NUM_ELEM_PER_UNIT ), (foundAllSize / nUnits) + 1 );
    
    
    #ifdef DEBUG
      sleep_for(std::chrono::milliseconds(SLEEP_TIME));
      cout << "ideal number of elements per unit: " << fmt( ideal, FRED ) << endl;
      
      cout << "Histogram: ";
      for( size_t i = 0; i < histo.lsize(); ++i ) {
        if( histo.local[i] ) cout << fmt( i, FCYN ) << ":" << fmt( histo.local[i], FRED ) << ", ";
      } 
      cout << endl;
    #endif
    

    vector<unitRange>::iterator  uRIt;
    
    // initialze with 0 (not sure if uchars are initialized with 0 by default, that's why i am doing this)
    for( uRIt = distr.begin( ); uRIt < distr.end( ); ++uRIt ) {
      uRIt->begin = 0;
      uRIt->end   = 0;
    }
    
    uRIt = distr.begin();
    uint acc = 0;
    size_t i;
    
    // actual calculation of distribution
    for( i = 0; i < histo.lsize(); ++i ) {
      acc += histo.local[i];
      
      if( acc >= ideal ){       
      
        uRIt->end = i;
        
        if( i+1 <= MAX_KEY ) (++uRIt)->begin = i+1;
        acc = 0;
      }
    }
    if( uRIt->begin > 0 ) uRIt->end = MAX_KEY;
  } // end of unit 0 only part
  
  // convert own team unit ID into team unit ID with 0
  dash::team_unit_t TeamUnit0ID = team.myid( );
  TeamUnit0ID.id = 0;

  dart_ret_t ret = dart_bcast(
                      static_cast<void*>( distr.data( ) ),  // buf 
                      distr.size( ) * sizeof(unitRange)  ,  // nelem
                      DART_TYPE_BYTE                     ,  // dtype
                      TeamUnit0ID                        ,  // root
                      team.dart_id( )                       // team
                   );
  
  if( DART_OK != ret ) cout << "An error while BCAST has occured!" << endl; 

  #ifdef DEBUG
    sleep_for(std::chrono::milliseconds(SLEEP_TIME));
    cout << "#" << fmt( myid, FBLUE, 2 ) << ": ";
    for( auto i : distr ) {
        cout << "Range: " 
          << fmt( i.begin, FYEL ) 
          << "-" 
          << fmt( i.end, FGREEN )
          << ", ";
      }
      cout << endl;
  #endif
}



int main( int argc, char* argv[] )
{  
  dash::init( &argc,&argv );
  int myid = static_cast<int>( dash::Team::GlobalUnitID( ).id );
  
  
  
  uint nrows = static_cast<uint>( atoi( argv[1] ) );
  uint ncols = static_cast<uint>( atoi( argv[2] ) );
  uint nelts = static_cast<uint>( atoi( argv[3] ) );
  
  winnow( nrows, ncols, nelts, myid );

  dash::finalize( );
}


