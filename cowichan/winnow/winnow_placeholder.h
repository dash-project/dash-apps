template< typename T, typename X, typename Y>
void winnow(
  uint const   nrows      ,
  uint const   ncols      ,
     T const & randMat    ,
     X const & threshMask , 
  uint const   nelts      ,
     Y       & points     )
{
  vector<pair<int, pair<int, int> > > values;
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
    Point pTmp = points[i];
    int index = i * chunk;
    pTmp.row = values[index].second.first;
    pTmp.col = values[index].second.second;
    points[i] = pTmp;
  }

}

