#include <vector>
#include <algorithm>

using std::vector;
using std::pair;
using std::make_pair;

using Point = struct{ uint row, col;};
using uint  = unsigned int ;


template< typename T, typename X>
void winnow(
                        uint const   nrows      ,
                        uint const   ncols      ,
                           T const & randMat    ,
                           X const & threshMask , 
                        uint const   nelts      ,
  vector<pair<POI_T, POI_T>>       & points     ){
  vector< pair<MATRIX_T, pair<POI_T, POI_T>> > values;
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      if (threshMask[i][j]) {
        values.push_back(make_pair(randMat[i][j], make_pair(i, j)));
      }
    }
  }
  sort(values.begin(), values.end());

  size_t n = values.size();
  size_t chunk = n / nelts;

  for (int i = 0; i < nelts; i++) {
    int index = i * chunk;
    points[i] = values[index].second;
  }

}

