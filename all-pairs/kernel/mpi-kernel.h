#ifndef MPI_KERNEL_H
#define MPI_KERNEL_H

#include <mpi.h>
#include "all-pairs-kernel.h"

class MPIKernel : public AllPairsKernel {

protected:
int      local_slots;
int    * send_data;
int    * recv_data;
size_t   sr_size;
int      repeat = 0;

MPIKernel(int internal_repeats = 1, std::string name = "MPI")
	: AllPairsKernel(internal_repeats, name)
	{}

~MPIKernel(){
	dash::barrier();
	MPI_Free_mem(send_data);
	MPI_Free_mem(recv_data);
}

void init(int repeats){
  sr_size = sizeof(int) * repeats * int_repeats;
  MPI_Alloc_mem(sr_size, MPI_INFO_NULL, &send_data);
  MPI_Alloc_mem(sr_size, MPI_INFO_NULL, &recv_data);

	MPI_Barrier(MPI_COMM_WORLD);	
}

public:
void reset(){
  repeat = 0;
}

};
#endif // MPI_KERNEL_H
