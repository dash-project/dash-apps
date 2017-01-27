#include <libdash.h>

#define DEBUG
#define SLEEP_TIME 30
#define MAX_KEY    100

#include <chrono>
#include <thread>

using std::cout;
using std::endl;
using std::vector;

using uint       = unsigned  int;
using uchar      = unsigned char;

using MTRX_TYPE  = uchar;

using Point      = struct{ MTRX_TYPE value; uint row, col;};
using pointRange = struct{ Point * curr, * end;           };


template< typename T = MTRX_TYPE >
inline void winnow(uint nrows, uint ncols, uint nelts, int myid);


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
    return os << "(" << static_cast<uint>( p.value ) << "," << p.row << "," << p.col << ")-";
  #else
    return os << p.row << " " << p.col << endl;
  #endif
}

int main( int argc, char* argv[] )
{  
  dash::init( &argc,&argv );
  int myid = static_cast<int>( dash::Team::GlobalUnitID( ).id );
  
  
  if( argc != 4 ){
    if( 0 == myid ){ cout << "3 Parameters expected!"              << endl
                          << "Usage: cowichan_winnow nRows nCols nElements" << endl
                          << "Then enter the matrix and the mask." << endl;
    }
    dash::finalize( );
    return 0;
  }
  
  uint nrows = static_cast<uint>( atoi( argv[1] ) );
  uint ncols = static_cast<uint>( atoi( argv[2] ) );
  uint nelts = static_cast<uint>( atoi( argv[3] ) );
  
  winnow( nrows, ncols, nelts, myid );

  dash::finalize( );
}


template<typename T = MTRX_TYPE>
inline void winnow(uint nrows, uint ncols, uint nelts, int myid){

  dash::NArray<T   , 2> matrix( nrows, ncols );
  dash::NArray<bool, 2> mask  ( nrows, ncols );
  
  dash::Team & team = dash::Team::All();
  size_t nUnits     = team.size();

  dash::Array<size_t> found( nUnits );
  vector<uint> histo( MAX_KEY );
  
  
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
  
  vector<Point> pointsLocal;
  // extra scope for loop (needed because of multiple declarations)
  {
    #ifdef DEBUG  // print error message if mask's and matrix's local size aren't identical
      if( mask.lend( ) - mask.lbegin( ) != matrix.lend( ) - matrix.lbegin( )) {
        cout << "On unit " << myid << " the local sizes of matrix and mask differ!\naborted on this unit" << endl;
        return;
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
      }
    }
    
    found[myid] = pointsLocal.size();    
    
    #ifdef DEBUG  // print pointslocal found size
      cout << "#" << myid << ": found :" << pointsLocal.size() << endl;
    #endif
  }
  
  
  sort( pointsLocal.begin( ), pointsLocal.end( ) );
  
  #ifdef DEBUG  // print points found local after sort
    std::this_thread::sleep_for(std::chrono::milliseconds(SLEEP_TIME)); 
    cout << "#" << myid << ": ";
    for( auto it : pointsLocal) {
      cout << it;
    } cout << endl;
  #endif

  
  vector<size_t> foundLclCpy( nUnits );  
  dash::barrier( );
  dash::copy(found.begin( ), found.end( ), foundLclCpy.data( ) );
  
  
  size_t foundAllSize = 0;
  for( size_t * i = foundLclCpy.data( ); i < foundLclCpy.data( ) + nUnits; ++i ) {
    foundAllSize += *i;
  }
  #ifdef DEBUG  // print foundAllSize
    std::this_thread::sleep_for(std::chrono::milliseconds(SLEEP_TIME)); 
    if( 0 == myid ) cout << "foundAllSize:" << foundAllSize << endl;
  #endif
  
  // will hold all points of all units in the end
  vector<Point> recvbuf( foundAllSize );
  vector<size_t> offset( nUnits + 1   );
  constexpr size_t pSize = sizeof( Point );
  
   {
    // initialze variables for loop
    vector<size_t> recvBytesPerUnit( nUnits );
    vector<size_t>::iterator fndIt  = foundLclCpy.begin( );
    vector<size_t>::iterator recvIt = recvBytesPerUnit.begin( );
    vector<size_t>::iterator offIt  = offset.begin          ( );
    
    
    /* initialze recvBytesPerUnit by calculating how much bytes are to be received from each unit
      * and initialze receive displacement/offset buffer 
      */
    *offIt = 0;
    ++offIt;
    
    for( ; recvIt != recvBytesPerUnit.end( ); ++recvIt, ++fndIt, ++offIt ) {
      *recvIt = *fndIt * pSize;
      *offIt  = *(offIt - 1) + *recvIt;
    }
    
    
    #ifdef DEBUG  // print recvbuf and offset
      std::this_thread::sleep_for(std::chrono::milliseconds(SLEEP_TIME)); 
      cout << "#" << myid << " recvbuf: ";
      for( auto it : recvBytesPerUnit ) {
        cout << it << ",";
      } cout << endl;
      
      
      std::this_thread::sleep_for(std::chrono::milliseconds(SLEEP_TIME)); 
      cout << "#" << myid << " offset: ";
      for( auto it : offset ) {
        cout << it << ",";
      } cout << endl;
      
    #endif
     // gather all data on every unit
     dart_ret_t ret = dart_allgatherv
                    (
                      static_cast<void*>( pointsLocal.data( ) ), // sendbuf
                      pointsLocal.size( ) * pSize              , // nsendelem
                      DART_TYPE_BYTE                           , // dtype
                      static_cast<void*>( recvbuf.data    ( ) ), // recvbuf
                      recvBytesPerUnit.data               ( )  , // nrecvelem
                      offset.data                         ( )  , // recvdispls
                      team.dart_id                        ( )    // teamid 
                    );
    
    if( DART_OK != ret ) cout << "An error has occured!" << endl; 
  
  }

  // unit 0 reorders the elements
  if( 0 == myid ) {
    
    nelts = foundAllSize;
    
    #ifdef DEBUG  // print all received elements
    std::this_thread::sleep_for(std::chrono::milliseconds(SLEEP_TIME));
    cout << "received data:" << endl
    
    for (int i = 0; i < foundAllSize; i++) {
      cout << recvbuf[i];
    } cout << endl;
    #endif
    
    // create final "array" and the pointer array (poinRanges) for reordering the recvbuf elements
    vector<Point> reordered( foundAllSize );
    vector<pointRange> pointRanges( nUnits );
    
    // initialze loop variables
    
    
    void *   rcvData         = recvbuf.data     ( );
    size_t * offPtr          =  offset.data ( ) + 1;
    pointRange * poiRangePtr = pointRanges.data ( );
    
    poiRangePtr->curr = static_cast<Point*>( rcvData );
    ++poiRangePtr;
    
    
    for( ; poiRangePtr < pointRanges.data ( ) + nUnits ; ++offPtr, ++poiRangePtr) {
      poiRangePtr->curr      = static_cast<Point*>( rcvData + *offPtr );
      (poiRangePtr - 1)->end = poiRangePtr->curr;
    } (poiRangePtr - 1)->end = static_cast<Point*>( rcvData + *offPtr);
    
    cout << endl << "here it comes" << endl;
    
    /* for( auto i : pointRanges ) {
      cout << "curr:" << *i.curr << ", end:" << *i.end << endl;
    } */
    
    
    
   /*  for(  ) {
      i.curr = & recvbuf[2];
      cout << *i.curr << endl;
    } */
    
    
  }
  
}





















