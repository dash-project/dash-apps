#ifndef RMA_KERNEL_H 
#define RMA_KERNEL_H 

#include <string>
#include <mpi.h>
#include "mpi-kernel.h"

/**
 * AllPairs Kernel concept
 */
class RMAKernel : public MPIKernel {

protected:

MPI_Win  window_send;
MPI_Win  window_recv;

public:

RMAKernel(int internal_repeats = 1, std::string name = "RMA")
  : MPIKernel(internal_repeats, name)
  {}

~RMAKernel(){
    dash::barrier();
    MPI_Win_unlock_all(window_recv);
    MPI_Win_unlock_all(window_send);  
    MPI_Win_free(&window_recv);
    MPI_Win_free(&window_send);
}

void init(int repeats){
		MPIKernel::init(repeats);
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
}

#endif // RMA_KERNEL_H 
