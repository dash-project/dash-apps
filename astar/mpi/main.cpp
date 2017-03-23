#include <mpi.h>
#include "JUtils.h"
#include "Astar.h"

int main (int argc, char* argv[]) {
	MPI_Init(NULL, NULL);
	
	JDurationManager dm;
	{
		Astar a;
		dm.start();
		
		a.run(1000);
		
		dm.stop();
    a.print_all(true);
    
    if (a.get_rank() == 0) {
      dm.print();
    }
	}
	MPI_Finalize();	
	return 0;
}
