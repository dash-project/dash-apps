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
   /* "gRow" represents the global row number of the local matrix data
    * the first local row has the initial global row number of "gRow"
    * "end" holds the global row number exakt one past the last row
    * number which is local at this unit.
    * "matP" will be used to linear iterate over the local data
    * "matBegin" will be used to access local data via []operator
    */
   auto gRow =        matOut.pattern().global({0,0})[0];
   auto end  = gRow + matOut.pattern().local_extents()[0];

   double nmax;
   value zero = {0,0};

   auto matBegin = matOut.lbegin();
   auto matP = matOut.lbegin();


   for( decltype(gRow) i = 0; gRow < end; ++gRow, ++i ) {

    nmax = 0;

    for( decltype(gRow) j = 0; j < nelts; ++j,++matP ) {
      if( gRow != j) {
        *matP = distance(points[gRow], points[j]);
        nmax  = max( nmax, *matP );
      }
    }

    matBegin[i*nelts+gRow] = nelts * nmax;
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