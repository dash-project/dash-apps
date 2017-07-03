// #ifndef DASH_ENABLE_LOGGING
// #define DASH_ENABLE_LOGGING
// #endif

#include <libdash.h>
#include <chrono>
#include <thread>
#include "../Terminal_Color.h"

#define DEBUG

#define MAX_KEY               99
#define MIN_NUM_ELEM_PER_UNIT 10


using std::cout   ;
using std::cin    ;
using std::endl   ;
using std::vector ;
using std::pair   ;

using dash::Team   ;
using dash::Array  ;
using dash::NArray ;
using dash::Shared ;
using dash::fill   ;

using uint       = unsigned  int ;
using uchar      = unsigned char ;
using POI_T      =          int  ;  //this type musst be signed!
using MATRIX_T   =         uchar ;

using Point      = struct{ MATRIX_T value; uint row, col; };
using pointRange = struct{ Point * curr, * end;           };  // unused until now
using unitRange  = struct{ MATRIX_T begin, end;           };  // unitValueRange


// typedef dash::CSRPattern< 1, dash::ROW_MAJOR, int > pattern_t;
// typedef typename pattern_t::size_type extent_t;

// static variables
static struct InputPar { uint nrows, ncols; } in;
static uint   nelts;
static int    myid ;

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
    return os << p.row << " " << p.col << "\n";
  #endif
}

#ifdef DEBUG
  using std::this_thread::sleep_for;
  inline void __sleep( uint const baseDur = 0, uint const mult = 10 )
  {
     uint SLEEP_TIME__ = (myid + 1) * mult + baseDur;
      sleep_for(std::chrono::milliseconds(SLEEP_TIME__)); 
  }
#endif

/*
 * One unit has the job to read in the parameters.
 * Because there's always a unit0, it reads the input parameter and
 * distributes them to the rest of the units.
 */
inline void ReadRowsNCols( )
{
  Shared<InputPar> input_transfer;

  if(0 == myid)
  {
    cin >> in.nrows;
    cin >> in.ncols;

    input_transfer.set(in);
  }
  input_transfer.barrier();
  in = input_transfer.get();
}


template< typename T = MATRIX_T >
inline void ReadMatricesAndNelts( NArray<T,2>& randMat, NArray<bool,2>& threshMask )
{
  Shared<uint> nelts_transfer;

  if(0 == myid)
  {
    //read matrices
    T tmp;
      for ( auto i : randMat ){
        scanf( "%u", &tmp )  , i = tmp;
      }
      bool tmpB;
      for ( auto i : threshMask ){
        scanf( "%u", &tmpB ) , i = tmpB;
      }
    cin >> nelts;

    nelts_transfer.set(nelts);
  }
  nelts_transfer.barrier();
  nelts = nelts_transfer.get();
}


template<typename T = MATRIX_T, typename X = pair<POI_T, POI_T> >
inline void winnow(
              uint const   nrows      ,
              uint const   ncols      ,
    NArray<   T,2> const & randMat    ,
    NArray<bool,2> const & threshMask ,
              uint const   nelts      ,
         vector<X>       & result     )
{  
  Team & team   = dash::Team::All ( );
  size_t nUnits = team.size       ( );
  
  if( 0 == myid ) cout << "nUnits: " << nUnits << endl;
  
  /* create global histo array for sorting
   * size += 1 for direct Index access -> histo[2]++ counts for value 2
   * size += 1 for additional value at the end of the 
   * histogram -> used for counter how many values were found
   */
  Array<uint> histo( (MAX_KEY + 2) * nUnits, dash::BLOCKED );
  // fill( histo.begin( ), histo.end( ), 0 );
  for(uint * it = histo.lbegin( ); it < histo.lend( ); ++it){*it = 0;}
  
  // local found points are gathered in this vector
  vector<Point> pointsLocal;
    
  // returns a object with the global row and column of the the local coordinates {0,0}
  auto globIndex = randMat.pattern( ).global( {0,0} );
  
  uint         gRow   = globIndex[0]           ;
  uint         gCol   = globIndex[1]           ;
  uint       & found  = histo.local[MAX_KEY+1] ;
     T const * matrEl = randMat.lbegin( )      ;
  
  
 /* read in local part of mask - matrix combination
  * and while doing that:
  * generate histogram and count the values (-> found)
  */
  for ( bool const * maskEl = threshMask.lbegin( );  maskEl < threshMask.lend( );  ++maskEl , ++matrEl, ++gCol )
  {
    if( gCol == ncols ) gCol = 0, ++gRow; // end of row -> next row in matrix/mask
    if( *maskEl       )
    {
      pointsLocal.push_back( Point{ *matrEl, gRow, gCol } );
      ++histo.local[*matrEl];
    }
  }
  found = pointsLocal.size();
  
  #ifdef DEBUG  // print points found local
    dash::barrier();  // only needed for better IO Output
    __sleep( );
    
    cout << "#" << fmt( myid, FBLUE, 2 ) << ": ";
    for( auto it : pointsLocal){ cout <<   it; }
    cout << endl;
    
    __sleep(20);
    cout << "#" << fmt( myid, FBLUE, 2 ) << ": found via pointsLocal:" << fmt( pointsLocal.size()   , FRED, 2 )  << "\n";
    cout << "#" << fmt( myid, FBLUE, 2 ) << ": found via histogram  :" << fmt( found                , FRED, 2 )  << "\n";
    // cout << "#" << fmt( myid, FBLUE, 2 ) << ": histogram.lsize:"       << fmt(found - histo.lbegin(), FRED, 2 )  << endl;
  #endif

  // is here a barrier needed? -> to wait for completion on unit0?!
  if( 0 != myid ) {
    dash::transform<uint>(
      histo.lbegin     ( ) ,
      histo.lend       ( ) ,
      histo.begin      ( ) ,
      histo.begin      ( ) ,
      dash::plus<uint> ( ) );
  }  
  
  
  // in this array will be the distribution info for creating buckets
  vector<unitRange> distr( nUnits );


  // unit 0 have to wait for rest to finish adding their values
  dash::barrier();  
  
  if( 0 == myid ) { // calculate bucket distribution

    /* calculate how much elements each unit should hold ideally
     * increment by one for safety garuantees in distribution
     * that is the assumption (ideal * nUnits > foundAllSize) == true
     */
     uint ideal = std::max( static_cast<size_t>( MIN_NUM_ELEM_PER_UNIT ), (found / nUnits) + 1 );
     
    #ifdef DEBUG
       __sleep();
      cout << "ideal number of elements per unit: " << fmt( ideal, FRED ) << endl;
      
      cout << "Histogram: ";
      for( size_t i = 0; i < histo.lsize(); ++i ) {
        if( histo.local[i] ) cout << fmt( i, FCYN ) << ":" << fmt( histo.local[i], FRED ) << ", ";
      } 
      cout << endl;
    #endif
     
  
    vector<unitRange>::iterator  uRIt;
    
   /* the loop for calculation of the distribution got a bit more complex
    * because i wanted to iterate only once over the "distr" vector
    * therefore no initialization beforehand is needed
    */ 
    // begin - 1 for loop logic (starting with ++)
    uRIt = distr.begin() - 1; 
    uint   acc    = 0;
       T   begin  = 0;
    uint * hisIt  = histo.lbegin();
    
    // actual calculation of distribution
    for( size_t i = 0; i < histo.lsize()-1 ; ++i ) {
      acc += *hisIt++;
      
      if( acc >= ideal ){
        
        (++uRIt)->begin = begin;
        uRIt->end = i;
        
        begin = i+1;
        acc   = 0;
      }
    }
    uRIt->end = MAX_KEY;
    
    // set the rest to zero
    while( ++uRIt < distr.end() ) { uRIt->begin = 0; uRIt->end = 0; }
    
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
    __sleep();
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
  
  
  
  
  /* each unit creates buckets that are send to the responsible unit later */
  //createBuckets( myid,team )
  
  
}



int main( int argc, char* argv[] )
{  
  dash::init( &argc,&argv );
  
  myid = dash::myid( );

  ReadRowsNCols( );
  
  NArray< MATRIX_T, 2 > randMat    ( in.nrows, in.ncols );
  NArray< bool    , 2 > threshMask ( in.nrows, in.ncols );

  #ifdef DEBUG  // print error message if mask's and matrix's local size aren't identical
    if( threshMask.local_size() != randMat.local_size() )
    {
      cout << "On unit " << myid
           << " the local sizes of matrix and mask differ!\naborted on this unit\n";
      return -1;
    }
  #endif
  
  ReadMatricesAndNelts( randMat, threshMask );
  
  vector< pair<POI_T, POI_T> > result_points(nelts);
  
  winnow( in.nrows, in.ncols, randMat, threshMask, nelts, result_points );

  dash::finalize( );
}


