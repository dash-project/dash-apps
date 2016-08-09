#ifndef MPI_ASYNC_KERNEL_H 
#define MPI_ASYNC_KERNEL_H 

#include <string>
#include <mpi.h>
#include "rma-kernel.h"

/**
 * AllPairs RMA Get Kernel 
 */
class MPIASyncKernel : public MPIKernel {

private:
MPI_Request    request_send;                                                  
MPI_Request    request_recv;   

public:

MPIASyncKernel(int internal_repeats = 1)
  : MPIKernel(internal_repeats, "MPI_ASYNC")
  {}

/**
 Perform repeated measures on given data point
 */
void run(int send, int recv){
      for(int r=0; r<int_repeats; r++){
        int sr_addr = repeat * int_repeats + r; // send / recieve addr
        if(myid == send){
          MPI_Isend (&(send_data[sr_addr]), 1, MPI_INT, recv, 99,    
                MPI_COMM_WORLD, &request_send);                                 
        }
        if(myid == recv){
          MPI_Irecv (&(recv_data[sr_addr]), 1, MPI_INT, recv, 99,    
                MPI_COMM_WORLD, &request_recv);
        }
        // Wait for request to complete                                   
        MPI_Wait (&request_recv, MPI_STATUS_IGNORE);
      }
      ++repeat;
}

};

#endif // MPI_ASYNC_KERNEL_H 
