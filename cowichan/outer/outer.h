inline double sqr(double const x) { return x * x; }


//calculates the distance between two points.
template<typename T = POI_T>
inline double distance(const pair<T, T>& x, const pair<T, T>& y) {
  return sqrt(sqr(x.first - y.first) + sqr(x.second - y.second));
}


template<typename T, typename X, typename Y>
inline void Outer(
      T const & points,
      X       & matOut,
      Y       & vec
    ){
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
    vec.local[i] = distance(make_pair(0,0), points[gRow]);
  }
  
  dash::barrier( );
}

