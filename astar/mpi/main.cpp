#include <mpi.h>
//#include <stdio.h>
//#include <vector>
//#include <map>
//#include <algorithm>
//#include "Puzzle.h"
#include "JUtils.h"
#include "Astar.h"
//#include <unistd.h>
//#include <cstddef>

int main (int argc, char* argv[]) {
	MPI_Init(NULL, NULL);
	
	JDurationManager dm;
	{
		dm.start();
		Astar a;
		a.run();
		
		dm.stop();
		dm.print();
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Finalize();	
	return 0;
}
