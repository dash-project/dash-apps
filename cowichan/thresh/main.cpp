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

  dash::NArray<T, 2> matSrc(nrows, ncols);
  //dash::NArray<T, 2> mask  (nrows, ncols);  //not used yet
  
  //read Input Matrix from stdin
  T tmp;
  if(0 == myid){
    for(auto i : matSrc){
      scanf("%d", &tmp);
      i = tmp;
    }
  }
  dash::barrier();
  
  //find max value in matSrc
  auto maxGlobIt = dash::max_element(matSrc.begin(), matSrc.end());
  T max          = *maxGlobIt;
  T maxPO        = max + 1; //max plus one
  
  //get number of Units running
  size_t num_units = dash::size();
  
  //calculate how many numbers are hold of each unit -> remark: each unit hold the number of every unit (rounding up)
  T histNumPerUnit = ((maxPO % num_units) > 0) ? (maxPO / num_units) + 1 : maxPO / num_units;
  
  //calculate tiling -> space on one unit to fit all numbers it's holding (last unit may be underfilled)
  size_t tiling    = histNumPerUnit * num_units;
  
  //for debugging
  size_t reqMem    = tiling * num_units;
 
  
  dash::Array<T> histo(reqMem);  //dash::TILE(tiling) -> not required anymore?
  
  //for debugging
/*   cout << "myid: " << myid << " has histo.lsize(): " << histo.lsize() << endl;
  
  for(uint i = 0; i < histo.lsize(); ++i){
    histo.local[i] = myid;
  }
  
  dash::barrier();
  
  if(0==myid){
    cout << endl << "--------" << endl;
    
    for(auto i : histo){
      cout << static_cast<uint>(i) << " ";
    }
    cout << endl << "--------" << endl;
  } */
  
  //zähle lokal die Vorkommen jeder Zahl und schreibe sie asynchron ins globale array
  for (T * i = matSrc.lbegin(); i < matSrc.lend(); ++i){
    if(*i){ //don't do anything if zero -> right now zeros aren't count at all
      int globIdx = (*i) * num_units + myid;
      histo.async[globIdx]++;      //should i make a differentiation between local and global?
    }
  }
  histo.async.push();
  
  
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
  
  dash::barrier();
  
  if(0 == myid) cout << endl 
                << "found max: "      << static_cast<uint>(max            ) << endl
                << "histNumPerUnit: " << static_cast<uint>(histNumPerUnit ) << endl
                << "tiling        : " << static_cast<uint>(tiling         ) << endl
                << "reqMem        : " << static_cast<uint>(reqMem         ) << endl;
  //cout << "myid: " << myid << " found max: " << static_cast<uint>(max) << endl;
  if(0 == myid) print2d(matSrc);
}





















