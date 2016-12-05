#include <libdash.h>
#include <iostream>
#include <cstdio>

using std::cout;
using std::endl;

typedef unsigned int  uint;
typedef unsigned char uchar;


template<typename T>
void print2d(T& mat) {
  for( int i=0; i<mat.extent(0); i++ ) {
    for( int j=0; j<mat.extent(1); j++ ) {
      cout << std::setw(3) << static_cast<uint>(mat(i,j))<< " ";
    }
    cout << endl;
  }
}

template<typename T>
uchar thresh(uint nrows, uint ncols, uint percent, T& matSrc, T& mask){
  //find max value in matSrc
  auto maxGlobIt = dash::max_element(matSrc.begin(), matSrc.end());
  uchar max = *maxGlobIt;
  
  uchar* histogram = new uchar[max+1];
  
  /* TO DO...
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      histogram[matrix[i][j]]++;
    }
  }

  int count = (nrows * ncols * percent) / 100;

  int prefixsum = 0;
  int threshold = nmax;

  for (int i = nmax; i >= 0 && prefixsum <= count; i--) {
    prefixsum += histogram[i];
    threshold = i;
  }

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      mask[i][j] = matrix[i][j] >= threshold;
    }
  }*/
  
  delete[] histogram;
  return max;
}

int main(int argc, char* argv[])
{
  dash::init(&argc,&argv);
  
  auto myid = dash::myid();
  
  if(argc != 4){
    if (0==myid) cout << "3 Parameters expected!" << endl << "Usage:cowichan_thresh nrows ncols percentage" << endl;
    dash::finalize();
    return 0;
  }
  
  uint nrows   = static_cast<uint>(atoi(argv[1]));
  uint ncols   = static_cast<uint>(atoi(argv[2]));
  uint percent = static_cast<uint>(atoi(argv[3]));
  
  dash::NArray<uchar, 2> matSrc(nrows, ncols);
  dash::NArray<uchar, 2> mask  (nrows, ncols);
  
  uchar tmp;
  if(0 == myid){
    for(auto i : matSrc){
      scanf("%d", &tmp);
      i = tmp;
    }
  }
  
  cout << static_cast<uint>(thresh(nrows, ncols, percent, matSrc, mask)) << endl;
  
  dash::barrier();
  if( myid==0 ) print2d(matSrc);

  dash::finalize();
}























