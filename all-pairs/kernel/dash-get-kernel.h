#ifndef DASH_GET_KERNEL_H
#define DASH_GET_KERNEL_H

class DashGetKernel : public AllPairsKernel {

typedef dash::TilePattern<1>               tpattern_t;
typedef dash::Array<int, long, tpattern_t> tarray_t;

private:

tarray_t testarray;
long     blocksize;
int      temp;
int      repeat = 0;

public:

DashGetKernel(int internal_repeats = 1)
  : AllPairsKernel(internal_repeats, "DASH_GET")
  {}

~DashGetKernel(){
  testarray.deallocate();
  dash::barrier();
}

void init(int repeats){
  blocksize = repeats * int_repeats;
  long size = dash::size() * blocksize;
  testarray.allocate(size, dash::BLOCKED);
}

void run(int send, int recv){
      for(int r=0; r<int_repeats; r++){
        long sr_addr = send * blocksize + repeat * int_repeats + r;
        if(myid == send){
          temp = testarray[sr_addr];
        }
     }
     ++repeat;
}
};
#endif  // DASH_GET_KERNEL_H
