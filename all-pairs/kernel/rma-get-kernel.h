#ifndef RMA_GET_KERNEL_H 
#define RMA_GET_KERNEL_H 

#include <string>
#include <mpi.h>
#include "all-pairs-kernel.h"

/**
 * AllPairs Kernel concept
 */
class RMAGetKernel : public AllPairsKernel {

private:

int      local_slots;
int    * send_data;
int    * recv_data;
size_t   sr_size;
int      repeat = 0;
MPI_Win  window_send;
MPI_Win  window_recv;


public:

RMAGetKernel(int internal_repeats = 1)
  : AllPairsKernel(internal_repeats, "RMA_GET")
  {}

~RMAGetKernel(){
    dash::barrier();
    MPI_Win_unlock_all(window_recv);
    MPI_Win_unlock_all(window_send);  
    MPI_Win_free(&window_recv);
    MPI_Win_free(&window_send);

    MPI_Free_mem(send_data);
    MPI_Free_mem(recv_data);
}

void init(int repeats){
  sr_size = sizeof(int) * repeats * int_repeats;
  MPI_Alloc_mem(sr_size, MPI_INFO_NULL, &send_data);
  MPI_Alloc_mem(sr_size, MPI_INFO_NULL, &recv_data);

    MPI_Win_create(send_data,
                  sr_size,
                  sizeof(int),
                  MPI_INFO_NULL,
                  MPI_COMM_WORLD,
                  &window_send);
    MPI_Win_create(recv_data,
                  sr_size,
                  sizeof(int),
                  MPI_INFO_NULL,
                  MPI_COMM_WORLD,
                  &window_recv);

    MPI_Win_fence(0, window_send);
    MPI_Win_fence(0, window_recv);
    MPI_Win_lock_all(0, window_send);
    MPI_Win_lock_all(0, window_recv);
    MPI_Barrier(MPI_COMM_WORLD);
}

void reset(){
  repeat = 0;
}

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
