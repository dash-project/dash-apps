#ifndef RMA_PUT_KERNEL_H
#define RMA_PUT_KERNEL_H

#include <string>
#include <mpi.h>
#include "rma-kernel.h"

/**
 * AllPairs RMA Put Kernel 
 */
class RMAPutKernel : public RMAKernel {

public:

    RMAPutKernel(int internal_repeats = 1)
        : RMAKernel(internal_repeats, "RMA_PUT")
    {}

    /**
     Perform repeated measures on given data point
     */
    void run(int send, int recv)
    {
        if(dash::myid() == send) {
            for(int r=0; r<int_repeats; r++) {
                int sr_addr = repeat * int_repeats + r; // send / recieve addr
                MPI_Put (&(send_data[sr_addr]), 1, MPI_INT, recv, sr_addr, 1,
                         MPI_INT, window_recv);
                MPI_Win_flush(recv, window_recv);
            }
            ++repeat;
        }
    }

};

#endif // RMA_PUT_KERNEL_H 
