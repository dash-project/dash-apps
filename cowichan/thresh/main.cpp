#include <libdash.h>
#include <iostream>
#include <cstdio>

using std::cout;
using std::endl;

typedef unsigned int  uint;
typedef unsigned char uchar;

#define USED_TYPE uchar


template<typename T = USED_TYPE>
inline void thresh(uint nrows, uint ncols, uint percent, int myid);


template<typename T>
void print2d(T& mat) {
  for( int i=0; i<mat.extent(0); i++ ) {
    for( int j=0; j<mat.extent(1); j++ ) {
      cout << std::setw(3) << static_cast<uint>(mat(i,j))<< " ";
    }
    cout << endl;
  }
}


int main(int argc, char* argv[])
{  
  dash::init(&argc,&argv);
  int myid = static_cast<int>( dash::myid() );
  
  if(argc != 4){
    if (0==myid){ cout << "3 Parameters expected!" 
                 << endl
                 << "Usage:cowichan_thresh nrows ncols percentage"
                 << endl;
    }
    dash::finalize();
    return 0;
  }
  
  uint nrows   = static_cast<uint>(atoi(argv[1]));
  uint ncols   = static_cast<uint>(atoi(argv[2]));
  uint percent = static_cast<uint>(atoi(argv[3]));
  
  thresh(nrows, ncols, percent, myid);

  dash::finalize();
}


template<typename T = USED_TYPE>
inline void thresh(uint nrows, uint ncols, uint percent, int myid){

  dash::NArray<T   , 2> matSrc(nrows, ncols);
  dash::NArray<bool, 2> mask  (nrows, ncols);  //not used yet
  
  // read Input Matrix from stdin
  if (0 == myid){
    T tmp;
    for ( auto i : matSrc ){
      scanf( "%d", &tmp );
      i = tmp;
    }
  }
  
  // wait for initialization of the matrix bevore calculating maximum
  dash::barrier();
  
  // find max value in matSrc
  auto maxGlobIt = dash::max_element(matSrc.begin(), matSrc.end());
  T max          = *maxGlobIt;
  T maxPO        =  max + 1  ;   //max plus one
  
  // get number of Units running
  size_t num_units = dash::size();
  
  dash::Array<uint> histo(maxPO * num_units, dash::BLOCKED);
  
  for (T * i = matSrc.lbegin(); i < matSrc.lend(); ++i) {
    ++histo.local[*i];
  }
  
  /* FOR REVIEW: i think the barrier is necessary because if unit 0 is still calculating
     while another unit sarts with dash::transform there could be a race condition */
  dash::barrier();
  
  // add the values of the local histogram to the histogram of unit 0
  if (0 != myid) {
    dash::transform<uint>( histo.lbegin     () ,
                           histo.lend       () ,
                           histo.begin      () ,
                           histo.begin      () ,
                           dash::plus<uint> ());
  }
  
  // wait for all units to finish adding
  dash::barrier();
  
/*   alternative
    if (0 != myid) {
    // Overwrite local histogram result with result histogram from unit 0:
    dash::copy(histo.begin(),           // Begin of block at unit 0
               histo.begin() + maxPO,   // End of block at unit 0
               histo.lbegin());
  } */
  dash::Shared<int> threshold;
  threshold.set(max);
  
  if (0 == myid) {
    
    #if 0
      for (uint j = 0; j < histo.lsize(); ++j) {
        if(histo.local[j]) cout << std::setw(3)   << j 
                << " counted: " << histo.local[j] << endl;
      }
    #endif
    //cout << "haeee" << endl;
    uint count = (nrows * ncols * percent) / 100;
    uint prefixsum = 0;
    int i;
    
    // for (uint * i = histo.lend(); i >= histo.lbegin() && prefixsum <= count; i--) {
    for (i = max; i >= 0 && prefixsum <= count; --i) {
      prefixsum += histo.local[i];
    }
    
    threshold.set (++i);
    cout << "original thresh: " << i << endl;
  }
  
  // wait for Unit 0 to finish
  //threshold.flush();
  //dash::barrier();
  threshold.barrier();

  {
  int threshLclCpy = threshold.get();
  
    T * src  = matSrc.lbegin( );
    bool * i = mask.lbegin(   );
    
        *i =   *src >= threshLclCpy;      
    while (i < mask.lend()) {
      *++i = *++src >= threshLclCpy;
    }
  #if 1
    cout << myid << " got threshold: " << threshLclCpy << endl;
  #endif
  }
  
  dash::barrier();
  if(0 == myid) {
    cout << "----------------------------" << endl;
    print2d(mask);
  }
}





















