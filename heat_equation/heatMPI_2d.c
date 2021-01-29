/*
 *   Heat conduction demo program
 *
 *  solves the heat equation on a 2D grid
 *
 *  March 2009
 *  Matthias.Lieber@tu-dresden.de
 *  Tobias.Hilbrich@zih.tu-dresden.de
 *
 *  Adapted: Jan 2013
 *
 *  Modified by Tobias Baumann (Summer 2016)
 *  tobias.baumann@mailbox.tu-dresden.de
 *
 *  MPI version.
 */

#include "heatMPI_2d.h"
#include "minimon.h"
#include <string.h>

/******************************************************
 * Allocate the heatGrid and initialize all variables.
 ******************************************************/
void heatAllocate(heatGrid *grid) {
    int xsize = grid->local_size_x;
    int ysize = grid->local_size_y;

    grid->theta    = (double**) malloc (sizeof(double*)*(xsize+2));
    grid->thetanew = (double**) malloc (sizeof(double*)*(xsize+2));
    grid->theta   [0] = (double*) malloc (sizeof(double)*(xsize+2)*(ysize+2));
    grid->thetanew[0] = (double*) malloc (sizeof(double)*(xsize+2)*(ysize+2));

    for (int i = 0; i < xsize+2; i++)
    {
        grid->theta   [i] = grid->theta   [0]+i*(ysize+2);
        grid->thetanew[i] = grid->thetanew[0]+i*(ysize+2);

        for (int j = 0; j < ysize+2; j++)
        {
            grid->theta   [i][j] = 0.0;
            grid->thetanew[i][j] = 0.0;
        }
    }

    grid->dx = 1.0;
    grid->dy = 1.0;
    grid->k = 1.0;
}

/******************************************************
 * Deallocate the heatGrid.
 ******************************************************/
void heatDeallocate(heatGrid *grid) {
    free (grid->theta[0]);
    free (grid->theta);
    grid->theta = NULL;
    free (grid->thetanew[0]);
    free (grid->thetanew);
    grid->thetanew = NULL;
}

/******************************************************
 * Initialize the grid ones.
 ******************************************************/
void heatInitialize(heatGrid* grid) {
    for (int x=1; x <= 100; x++) {
        for (int y=1; y <= 100; y++) {
            grid->theta[x][y] = 1.0;
            grid->thetanew[x][y] = 1.0;
        }
    }
}

/******************************************************
 * Calculate the total energy (sum of theta in all grid cells)
 ******************************************************/
void heatTotalEnergy(heatGrid* grid, double *energy) {
  *energy = 0.0;
  for (int x=1; x <= grid->local_size_x; x++)
      for (int y=1; y <= grid->local_size_y; y++)
          *energy += grid->theta[x][y];
}

/******************************************************
 * Function to setup MPI data.
 *
 * (1) Initializes MPI
 * (2) Creates a cartesian communicator for border exchange
 * (3) Tell each process where is is located an how many cells it has to handle
 * (4) Sets up helpful data-type and MPI buffer
 *
 ******************************************************/
void heatMPISetup (heatGrid *grid, dataMPI* configMPI) {
    int dims[2] = {0,0},
        periods[2] = {1,1},
        coords[2];
    int buf_size;
    char *buf;

    /* ==== (1) ==== */
    /* Base init*/
    MPI_Comm_rank (MPI_COMM_WORLD, &configMPI->rank);
    MPI_Comm_size (MPI_COMM_WORLD, &configMPI->size);

    /* ==== (2) ==== */
    /* Create cartesian communicator*/
    MPI_Dims_create (configMPI->size, 2, dims);
    MPI_Cart_create (MPI_COMM_WORLD, 2, dims, periods, 0, &configMPI->cart);

    configMPI->dim_x = dims[0];
    configMPI->dim_y = dims[1];

    /* Store neighbors in the grid */
    MPI_Cart_shift (configMPI->cart, 0, 1, &configMPI->left, &configMPI->right);
    MPI_Cart_shift (configMPI->cart, 1, 1, &configMPI->up,   &configMPI->down );

    /* ==== (3) ==== */
    /* Create partitioning of overall grid on processes */
    MPI_Cart_coords (configMPI->cart, configMPI->rank, 2, coords); /*My coordinate*/

    /* Save size of cartesian array */
    configMPI->num_cells_x = dims[0];
    configMPI->num_cells_y = dims[1];

    /* Save position of this process in cartesian array */
    configMPI->coord_x = coords[0];
    configMPI->coord_y = coords[1];

    /* Save how many cells a process has to handle */
    grid->local_size_x = grid->global_size_x/dims[0];
    grid->local_size_y = grid->global_size_y/dims[1];

    /* Recalc global size because of rounding */
    grid->global_size_x = grid->local_size_x*dims[0];
    grid->global_size_y = grid->local_size_y*dims[1];

    /* ==== (4) ==== */
    /* Create datatype to communicate one row */
    MPI_Type_vector (
            grid->local_size_x, /* #blocks */
            1, /* #elements per block */
            grid->local_size_y+2, /* #stride */
            MPI_DOUBLE, /* old type */
            &configMPI->rowtype /* new type */ );
    MPI_Type_commit (&configMPI->rowtype);
}

/******************************************************
 * Function to free and finalize MPI.
 ******************************************************/
void heatMPIFree (dataMPI* configMPI) {
    MPI_Type_free (&configMPI->rowtype);
    MPI_Comm_free (&configMPI->cart);
}

void calcInnerGrid(heatGrid* grid, dataMPI* mympi, double dt)
{
    double dtheta;

    /* Calculate time step for inner cells: read from theta, write new timestep to thetanew. Dont calc outer cells. */
    for(int x = 2; x < grid->local_size_x; x++) {
        for(int y = 2; y < grid->local_size_y; y++) {
            dtheta = ( grid->theta[x-1][y] + grid->theta[x+1][y] - 2*grid->theta[x][y]) / (grid->dx * grid->dx)
                   + ( grid->theta[x][y-1] + grid->theta[x][y+1] - 2*grid->theta[x][y]) / (grid->dy * grid->dy);
            grid->thetanew[x][y] = grid->theta[x][y] + grid->k * dtheta * dt;
        }
    }
}

// calc inner cells of the row at position row, dont calc corners because they depend on other ghost cells
void calcRow(heatGrid* grid, dataMPI* mympi, double dt, int row)
{
    double dtheta;

    for(int x = 2; x < grid->local_size_x; x++) {
        dtheta = (grid->theta[x-1][row] + grid->theta[x+1][row] - 2*grid->theta[x][row]) / (grid->dx * grid->dx)
               + (grid->theta[x][row-1] + grid->theta[x][row+1] - 2*grid->theta[x][row]) / (grid->dy * grid->dy);
        grid->thetanew[x][row] = grid->theta[x][row] + grid->k * dtheta * dt;
    }
}

// calc the column at position column between startpoint and endpoint (used to calc with and without corners)
void calcColumn(heatGrid* grid, dataMPI* mympi, double dt, int column, int start, int end)
{
    double dtheta;

    for(int y = start; y < end; y++) {
        dtheta = (grid->theta[column-1][y] + grid->theta[column+1][y] - 2*grid->theta[column][y]) / (grid->dx * grid->dx)
               + (grid->theta[column][y-1] + grid->theta[column][y+1] - 2*grid->theta[column][y]) / (grid->dy * grid->dy);
        grid->thetanew[column][y] = grid->theta[column][y] + grid->k * dtheta * dt;
    }
}

MiniMon minimon{};

/******************************************************
 * Main program and time stepping loop.
 ******************************************************/
int main (int argc, char** argv) {
    minimon.enter();

    heatGrid mygrid;
    double dt, energyLocal, energyInitial, energyFinal;
    int step, nsteps;
    dataMPI mympi;

    /* initialize MPI */
    MPI_Init(&argc, &argv);

    if(argc != 4) {
        perror("Usage: heatMPI <x-size> <y-size> <iterations>");
        return 1;
    }

    mygrid.global_size_x = strtol(argv[1], NULL, 10);
    mygrid.global_size_y = strtol(argv[2], NULL, 10);

    /* setup MPI */
    heatMPISetup(&mygrid, &mympi);

    /* create heatGrid and initialize variables */
    heatAllocate(&mygrid);
    nsteps = strtol(argv[3], NULL, 10);
    energyInitial = energyFinal = 0.0;
    dt = 0.05;

    if (mympi.rank == 0)
        heatInitialize(&mygrid);

    /* energy of initial grid */
    heatTotalEnergy(&mygrid, &energyLocal);
    MPI_Reduce(&energyLocal, &energyInitial, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (mympi.rank == 0) {
        printf("Running Jacobi on %d processors with (%d, %d) elements.\n", mympi.size, mympi.num_cells_x, mympi.num_cells_y);
        printf("Array Dimensions: %d %d\n", mygrid.global_size_x, mygrid.global_size_y);
        printf("Block Dimensions: %d %d\n", mygrid.local_size_x, mygrid.local_size_y);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /*************/
    /* Main loop */
    /*************/
    minimon.enter();
    for(step=0 ; step<nsteps ; step++) {
        minimon.enter();

        //start communication
        minimon.enter();
        MPI_Request send_requests[4];
        MPI_Request recv_requests[4];

        //Irecv before Isend for better performance

        /*Receive Right border column from right neighbor*/
        MPI_Irecv(&(mygrid.theta[mygrid.local_size_x+1][1]),
                mygrid.local_size_y, MPI_DOUBLE, mympi.right, step, mympi.cart, &recv_requests[0]);
        /*Send left column to left neighbor*/
        MPI_Isend(&(mygrid.theta[1][1]),
                mygrid.local_size_y, MPI_DOUBLE, mympi.left, step, mympi.cart, &send_requests[0]);

        /*Receive Left border column from left neighbor*/
        MPI_Irecv(&(mygrid.theta[0][1]),
                mygrid.local_size_y, MPI_DOUBLE, mympi.left, step, mympi.cart, &recv_requests[1]);
        /*Send right column to right neighbor*/
        MPI_Isend(&(mygrid.theta[mygrid.local_size_x][1]),
                mygrid.local_size_y, MPI_DOUBLE, mympi.right, step, mympi.cart, &send_requests[1]);

        /*Receive lower border row from bottom neighbor*/
        MPI_Irecv(&(mygrid.theta[1][mygrid.local_size_y+1]),
                1, mympi.rowtype, mympi.down, step, mympi.cart, &recv_requests[2]);
        /*Send upper row to top neighbor*/
        MPI_Isend(&(mygrid.theta[1][1]),
                1, mympi.rowtype, mympi.up, step, mympi.cart, &send_requests[2]);

        /*Receive upper border row from top neighbor*/
        MPI_Irecv(&(mygrid.theta[1][0]),
                1, mympi.rowtype, mympi.up, step, mympi.cart, &recv_requests[3]);
        /*Send lower row to bottom neighbor*/
        MPI_Isend(&(mygrid.theta[1][mygrid.local_size_y]),
                1, mympi.rowtype, mympi.down, step, mympi.cart, &send_requests[3]);

        minimon.leave("calc update async");
        //overlap communication with computation
        minimon.enter();
        calcInnerGrid(&mygrid, &mympi, dt);
        minimon.leave("calc inner");

        minimon.enter();
        // wait for all Recv to complete and save switch-case statement
        MPI_Waitall(4, recv_requests, MPI_STATUSES_IGNORE);
        minimon.leave("calc wait");

        /* calcColumn has more performance than calcRow because of memory access pattern,
           thus using calcColumn to calc corners */

        minimon.enter();
        // calc right border and corners
        calcColumn(&mygrid, &mympi, dt, mygrid.local_size_x, 1, mygrid.local_size_y+1);
        // calc left border and corners
        calcColumn(&mygrid, &mympi, dt, 1,                   1, mygrid.local_size_y+1);
        // calc lower border
        calcRow(&mygrid, &mympi, dt, mygrid.local_size_y);;
        // calc upper border
        calcRow(&mygrid, &mympi, dt, 1);
        minimon.leave("calc bound");

        //Sends have to be finished before swaping pointers
        MPI_Waitall(4, send_requests, MPI_STATUSES_IGNORE);

        //swap pointer
        double **tmp;
        tmp = mygrid.theta;
        mygrid.theta = mygrid.thetanew;
        mygrid.thetanew = tmp;
        minimon.leave("calc iter");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    minimon.leave("calc total");

    /* Energy of final grid */
    heatTotalEnergy(&mygrid, &energyLocal);
    MPI_Reduce(&energyLocal, &energyFinal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    /* Work only for master process*/
    if (mympi.rank == 0) {

        printf("== Energy Conservation Check ==\n");
        printf("initial Energy: %20.3f\n", energyInitial);
        printf("  final Energy: %20.3f\n", energyFinal);
        printf("    Difference: %20.3f\n", energyFinal - energyInitial);
        printf("== Time Statistics ==\n");
        printf("Completed %d iterations.\n", nsteps);
    }

    heatDeallocate(&mygrid);

    /* Finalize MPI*/
   minimon.leave("total");

  if(mympi.rank == 0)
    std::cout << "# unit_id;function_name;num_calls;avg_runtime;min_runtime;max_runtime"
              << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);
   for(auto i = 0; i < mympi.size; ++i) {
    if(i == mympi.rank)
      minimon.print(mympi.rank);
    MPI_Barrier(MPI_COMM_WORLD);
  }
    heatMPIFree (&mympi);
    MPI_Finalize();

    return 0;
}
