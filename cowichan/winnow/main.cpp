#include <libdash.h>

//#define DEBUG

//#include <chrono>
//#include <thread>

using std::cout;
using std::endl;
using std::vector;

using uint      = unsigned  int;
using uchar     = unsigned char;

using MTRX_TYPE = uchar;

using Point = struct{ MTRX_TYPE value;
                      uint row, col;};
                      

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
  return os << "(" << static_cast<uint>( p.value ) << "," << p.row << "," << p.col << ")";
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

  dash::NArray <T   , 2> matrix( nrows, ncols );
  dash::NArray <bool, 2> mask  ( nrows, ncols );
  
  dash::Team & team = dash::Team::All();
  size_t nUnits     = team.size();

  return;
  dash::Array  < uint  > found ( nUnits );
  // dash::Array  < Point > points( nelts        );
  
  
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
  
  vector< Point > pointsLocal;
  // extra scope for loop (needed because of multiple declarations)
  {
    #ifdef DEBUG
      if( mask.lend( ) - mask.lbegin( ) != matrix.lend( ) - matrix.lbegin( )) {
        cout << "On unit " << myid << " the local sizes of matrix and mask differ!\naborted on this unit" << endl;
        return;
      }
    #endif
    
    auto globIndex = matrix.pattern( ).global( {0,0} );
    uint gRow  = globIndex[0];
    uint gCol  = globIndex[1];
    T * matrEl = matrix.lbegin( );
    
    for ( bool * maskEl = mask.lbegin( );  maskEl < mask.lend( );  ++maskEl , ++matrEl, ++gCol ) {
      if( gCol == ncols ) gCol = 0, ++gRow;
      if( *maskEl ) pointsLocal.push_back( Point{ *matrEl, gRow, gCol } );
    }
    
    found[myid] = pointsLocal.size();
    #ifdef DEBUG
      cout << "#" << myid << ": found :" << pointsLocal.size() << endl;
    #endif
  }
  
  sort( pointsLocal.begin( ), pointsLocal.end( ) );
  
  
  std::vector<uint> foundLclCpy(nUnits);
  
  dash::Shared<dash::GlobPtr<Point>> destUnit0;
  dash::barrier( );
  
  dash::copy(found.begin( ), found.end( ), foundLclCpy.data( ) );
  
  uint foundAllSize = 0;
  for( uint * i = foundLclCpy.data( ); i < foundLclCpy.data( ) + nUnits; ++i ) {
    foundAllSize += *i;
  }
    

    
/*   if( 0 == myid ) {
    //cout << "foundAllSize:" << foundAllSize << endl;
    destUnit0.set( dash::memalloc<Point>( foundAllSize ) );
  } */
  
  //destUnit0.barrier( );
  //dash::GlobPtr<Point> destUnit0LP = destUnit0.get( );
  
  
  //dash::GlobIter<Point>(destUnit0LP);
  //destUnit0LP += 10;
  
  // das ist der Teil an dem ich gerade scheitere
  //dash::copy( pointsLocal.begin( ), pointsLocal.end( ), destUnit0LP);
  
  
  //dash::barrier( );
  if( 0 == myid ) {
    cout << "----------------------------" << endl;
    //print result
  }
}





















