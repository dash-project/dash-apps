#ifndef ALL_PAIRS_KERNEL_H
#define ALL_PAIRS_KERNEL_H

#include <string>

/**
 * AllPairs Kernel concept
 */
class AllPairsKernel {

private:

const int         int_repeats;
const std::string kernel_name = "Demo";

public:

AllPairsKernel(int internal_repeats = 10)
  : int_repeats(internal_repeats)
  {}

~AllPairsKernel(){}

/**
 Perform repeated measures on given data point
 */
void run(int send, int recv){
  if(dash::myid() == send){
    std::cout << "(" << send << "," << recv << ")" << std::endl;
  }
}

const std::string name(){
  return this->kernel_name;
}

};



#endif // ALL_PAIRS_KERNEL_H
