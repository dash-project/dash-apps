#include <mpi.h>
#include "JUtils.h"
#include "Astar.h"

int main (int argc, char* argv[]) {
	MPI_Init(NULL, NULL);
	
	JDurationManager dm;
	{
		Astar a;
		dm.start();
		
		a.run();
		
		dm.stop();
		
		if (a.get_rank() == 0) {
			//a.print();
			dm.print();
		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		if (a.get_rank() == 1) {
			//a.print();
			dm.print();
		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		if (a.get_rank() == 2) {
			//a.print();
			dm.print();
		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		if (a.get_rank() == 3) {
      //a.print();
			dm.print();
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Finalize();	
	return 0;
}
