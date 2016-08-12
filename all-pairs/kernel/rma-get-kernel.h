#ifndef RMA_GET_KERNEL_H 
#define RMA_GET_KERNEL_H 

#include <string>
#include <mpi.h>
#include "rma-kernel.h"

/**
 * AllPairs RMA Get Kernel 
 */
class RMAGetKernel : public RMAKernel {

public:

RMAGetKernel(int internal_repeats = 1)
  : RMAKernel(internal_repeats, "RMA_GET")
  {}

/**
 Perform repeated measures on given data point
 */
void run(int send, int recv){
  if(dash::myid() == send){
      for(int r=0; r<int_repeats; r++){
        int sr_addr = repeat * int_repeats + r; // send / recieve addr
        MPI_Get(&(recv_data[sr_addr]), 1, MPI_INT, recv, sr_addr, 1,
                  MPI_INT, window_send);
        MPI_Win_flush(recv, window_send);
      }
      ++repeat;
  }
}

};

#endif // RMA_GET_KERNEL_H 
