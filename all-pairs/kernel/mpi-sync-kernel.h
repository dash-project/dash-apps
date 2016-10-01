#ifndef MPI_SYNC_KERNEL_H 
#define MPI_SYNC_KERNEL_H 

#include <string>
#include <mpi.h>
#include "rma-kernel.h"

/**
 * AllPairs RMA Get Kernel 
 */
class MPISyncKernel : public MPIKernel {

public:

MPISyncKernel(int internal_repeats = 1)
  : MPIKernel(internal_repeats, "MPI_SYNC")
  {}

/**
 Perform repeated measures on given data point
 */
void run(int send, int recv){
      if(send == recv){
        // skip as this is not possible
        return;
      }
      for(int r=0; r<int_repeats; r++){
        int sr_addr = repeat * int_repeats + r; // send / recieve addr
        if(myid == send){
            //std::cout << "SENDER: send: " << send << " recv: " << recv << std::endl;
            MPI_Send (&(send_data[sr_addr]), 1, MPI_INT, recv, 99, MPI_COMM_WORLD);
        }
        if(myid ==recv){
            //std::cout << "RECEIVER: send: " << send << " recv: " << recv << std::endl;
            MPI_Recv (&(recv_data[sr_addr]), 1, MPI_INT, send, 99, MPI_COMM_WORLD,
                    &status);     
        }
      }
      ++repeat;
}

};

#endif // MPI_SYNC_KERNEL_H 
