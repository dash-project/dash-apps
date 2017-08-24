using  std::max;
using  std::pair;
using  std::vector;
using  std::make_pair;

using dash::Team;
using dash::Array;
using dash::NArray;
using dash::barrier;
using dash::team_unit_t;


inline double sqr(double const x) { return x * x; }


//calculates the distance between two points.
inline double distance( value const & x, value const & y ){
  return sqrt(sqr(x.row - y.row) + sqr(x.col - y.col));
}


inline void Outer(
  vector< value      > const & points,
  NArray< double, 2  >       & matOut,
  Array < double     >       & vec   ,
                  uint         nelts )
{
   // comment grow, end i and j!!!
   uint gRow =        matOut.pattern().global({0,0})[0];
   uint end  = gRow + matOut.pattern().local_extents()[0];
   double nmax;
   value zero = {0,0};
   // zero.row = 0;
   // zero.col = 0;

   for( uint i = 0; gRow < end; ++gRow, ++i ) {
    nmax = 0;
    for( uint j = 0; j < nelts; ++j ) {
      if( gRow != j) {
        matOut.local[i][j] = distance(points[gRow], points[j]);
        nmax = max( nmax, static_cast<double>( matOut.local[i][j] ) );
      }
    }
    matOut.local[i][gRow] = nelts * nmax;
    vec.local[i] = distance( zero, points[gRow] );
  }
  
  barrier( );
}


template<typename T>
inline void BroadcastPointsToUnits( vector<T> & points )
{
  dart_ret_t ret = dart_bcast(
                      static_cast<void*>( points.data() ),  // buf 
                      points.size( ) * sizeof(T)         ,  // nelts
                      DART_TYPE_BYTE                     ,  // dtype
                      dash::team_unit_t(0)               ,  // root
                      dash::Team::All( ).dart_id( )      ); // team
                      
  if( DART_OK != ret ) cout << "An error while BCAST has occured!" << endl; 
}