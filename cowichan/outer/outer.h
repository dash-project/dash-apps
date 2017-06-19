using std::pair;
using std::make_pair;
using std::vector;
using std::max;
using dash::NArray;
using dash::Array;
using dash::team_unit_t;
using dash::Team;
using dash::barrier;

inline double sqr(double const x) { return x * x; }


//calculates the distance between two points.
template<typename T = POI_T>
inline double distance( pair< T, T > const & x, pair< T, T > const & y ){
  return sqrt(sqr(x.first - y.first) + sqr(x.second - y.second));
}


template< typename T = POI_T>
inline void Outer(
  vector< pair<T, T> > const & points,
  NArray< double, 2  >       & matOut,
  Array < double     >       & vec   ,
                  uint         nelts )
{
   // comment grow, end i and j!!!
   uint gRow =        matOut.pattern().global({0,0})[0];
   uint end  = gRow + matOut.pattern().local_extents()[0];
   double nmax;

   for( uint i = 0; gRow < end; ++gRow, ++i ) {
    nmax = 0;
    for( uint j = 0; j < nelts; ++j ) {
      if( gRow != j) {
        matOut.local[i][j] = distance(points[gRow], points[j]);
        nmax = max( nmax, static_cast<double>( matOut.local[i][j] ) );
      }
    }
    matOut.local[i][gRow] = nelts * nmax;
    vec.local[i] = distance( make_pair(0,0), points[gRow] );
  }
  
  barrier( );
}


template<typename T = POI_T>
inline void BroadcastPointsToUnits( vector< pair<T,T> > & points )
{
  team_unit_t TeamUnit0ID = Team::All( ).myid( );
  TeamUnit0ID.id = 0;
  dart_ret_t ret = dart_bcast(
                      static_cast<void*>( points.data( ) ),  // buf 
                      points.size( ) * sizeof(pair<T,T>)  ,  // nelts
                      DART_TYPE_BYTE                      ,  // dtype
                      TeamUnit0ID                         ,  // root
                      Team::All().dart_id( )              ); // team
                      
  if( DART_OK != ret ) cout << "An error while BCAST has occured!" << endl; 
}