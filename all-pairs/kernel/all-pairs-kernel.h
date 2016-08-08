#ifndef ALL_PAIRS_KERNEL_H
#define ALL_PAIRS_KERNEL_H

#include <string>

/**
 * AllPairs Kernel concept
 */
class AllPairsKernel {

protected:

const int         int_repeats;
const int         myid  = 0;
std::string kernel_name = "Demo";

AllPairsKernel(int internal_repeats, std::string name)
  : int_repeats(internal_repeats),
    myid(dash::myid()),
    kernel_name(name)
  {}

public:

AllPairsKernel(int internal_repeats = 1)
  : int_repeats(internal_repeats)
  {}

~AllPairsKernel(){}

/**
 * Initialize kernel
 */
void init(int repeats){}

/**
 * Reset kernel for new run
 */
void reset(){}

/**
 Perform repeated measures on given data point
 */
void run(int send, int recv){
  if(dash::myid() == send){
    std::cout << "(" << send << "," << recv << ")" << std::endl;
  }
}

const std::string getName(){
  return this->kernel_name;
}

int getInternalRepeats(){
  return this->int_repeats;
}
};



#endif // ALL_PAIRS_KERNEL_H
