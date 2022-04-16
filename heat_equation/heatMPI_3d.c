
#include "heatMPI_3d.h"
#include "minimon.h"
#include <string.h>

/******************************************************
 * Allocate the heatGrid and initialize all variables.
 ******************************************************/
void heatAllocate(heatGrid *grid) {
    long xsize = grid->local_size_x + 2;
    long ysize = grid->local_size_y + 2;
    long zsize = grid->local_size_z + 2;

    long total_size = (xsize)*(ysize)*(zsize);

    grid->theta    = (double***) malloc (sizeof(double**) * xsize);
    grid->thetanew = (double***) malloc (sizeof(double**) * xsize);
    grid->theta   [0] = (double**) malloc (sizeof(double*)*(xsize)*(ysize));
    grid->thetanew[0] = (double**) malloc (sizeof(double*)*(xsize)*(ysize));
    grid->theta   [0][0] = (double*) malloc (sizeof(double) * total_size);
    grid->thetanew[0][0] = (double*) malloc (sizeof(double) * total_size);

    int off_x = ysize;
    int off_y = zsize;
    long offset_x = 0;
    long offset_y = 0;
    for (int x = 0; x < xsize; x++)
    {
      grid->theta   [x] = grid->theta[0] + offset_x;
      grid->thetanew[x] = grid->thetanew[0] + offset_x;

      for (int y = 0; y < ysize; y++) {
        grid->theta   [x][y] = grid->theta[0][0] + offset_y;
        grid->thetanew[x][y] = grid->thetanew[0][0] + offset_y;
        for (int z = 0; z < zsize; z++)
        {
            grid->theta   [x][y][z] = 0.0;
            grid->thetanew[x][y][z] = 0.0;
        }
        offset_y += off_y;
      }
      offset_x += off_x;
    }

    grid->dx = 1.0;
    grid->dy = 1.0;
    grid->dz = 1.0;
    grid->k = 1.0;
}

/******************************************************
 * Deallocate the heatGrid.
 ******************************************************/
void heatDeallocate(heatGrid *grid) {
    free (grid->theta[0][0]);
    free (grid->theta[0]);
    free (grid->theta);
    grid->theta = NULL;
    free (grid->thetanew[0][0]);
    free (grid->thetanew[0]);
    free (grid->thetanew);
    grid->thetanew = NULL;
}

/******************************************************
 * Initialize the grid ones.
 ******************************************************/
void heatInitialize(heatGrid* grid) {
  for (int x=1; x <= grid->local_size_x; x++) {
    for (int y=1; y <= grid->local_size_y; y++) {
      for (int z=1; z <= grid->local_size_z; z++) {
        if(x <= 100 && y <= 100 && z <= 100) {
          grid->theta[x][y][z] = 1.0;
          grid->thetanew[x][y][z] = 1.0;
        }
      }
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
        for (int z=1; z <= grid->local_size_z; z++)
          *energy += grid->theta[x][y][z];
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
    int dims[3] = {0,0,0},
        periods[3] = {1,1,1},
        coords[3];

    /* ==== (1) ==== */
    /* Base init*/
    MPI_Comm_rank (MPI_COMM_WORLD, &configMPI->rank);
    MPI_Comm_size (MPI_COMM_WORLD, &configMPI->size);

    /* ==== (2) ==== */
    /* Create cartesian communicator*/
    MPI_Dims_create (configMPI->size, 3, dims);
    MPI_Cart_create (MPI_COMM_WORLD, 3, dims, periods, 0, &configMPI->cart);

    configMPI->dim_x = dims[0];
    configMPI->dim_y = dims[1];
    configMPI->dim_z = dims[2];

    /* Store neighbors in the grid */
    MPI_Cart_shift (configMPI->cart, 0, 1, &configMPI->up,   &configMPI->down );
    MPI_Cart_shift (configMPI->cart, 1, 1, &configMPI->left, &configMPI->right);
    MPI_Cart_shift (configMPI->cart, 2, 1, &configMPI->front,   &configMPI->back );

    /* ==== (3) ==== */
    /* Create partitioning of overall grid on processes */
    MPI_Cart_coords (configMPI->cart, configMPI->rank, 3, coords); /*My coordinate*/

    /* Save size of cartesian array */
    configMPI->num_cells_x = dims[0];
    configMPI->num_cells_y = dims[1];
    configMPI->num_cells_z = dims[2];

    /* Save position of this process in cartesian array */
    configMPI->coord_x = coords[0];
    configMPI->coord_y = coords[1];
    configMPI->coord_z = coords[2];

    /* Save how many cells a process has to handle */
    grid->local_size_x = grid->global_size_x/dims[0];
    grid->local_size_y = grid->global_size_y/dims[1];
    grid->local_size_z = grid->global_size_z/dims[2];

    /* Recalc global size because of rounding */
    grid->global_size_x = grid->local_size_x*dims[0];
    grid->global_size_y = grid->local_size_y*dims[1];
    grid->global_size_z = grid->local_size_z*dims[2];

    /* ==== (4) ==== */
    /* Create datatype to communicate one row */
    MPI_Type_vector (
            grid->local_size_y, /* #blocks */
            grid->local_size_z, /* #elements per block */
            grid->local_size_z+2, /* #stride */
            MPI_DOUBLE, /* old type */
            &configMPI->tb_type /* new type */ );
    MPI_Type_commit (&configMPI->tb_type);

    MPI_Type_vector (
            grid->local_size_x, /* #blocks */
            grid->local_size_z, /* #elements per block */
            (grid->local_size_y+2) * (grid->local_size_z+2), /* #stride */
            MPI_DOUBLE, /* old type */
            &configMPI->lr_type /* new type */ );
    MPI_Type_commit (&configMPI->lr_type);

    int sizes[3] = {grid->local_size_x+2, grid->local_size_y+2, grid->local_size_z+2};
    int sizes_sub[3] = {grid->local_size_x, grid->local_size_y,1};
    int starts[3] = {0,0,0};
    MPI_Type_create_subarray (
            3,
            sizes,
            sizes_sub,
            starts,
            MPI_ORDER_C, /* #blocks */
            MPI_DOUBLE, /* old type */
            &configMPI->fb_type /* new type */ );
    MPI_Type_commit (&configMPI->fb_type);
}

/******************************************************
 * Function to free and finalize MPI.
 ******************************************************/
void heatMPIFree (dataMPI* configMPI) {
    MPI_Type_free (&configMPI->tb_type);
    MPI_Type_free (&configMPI->lr_type);
    MPI_Type_free (&configMPI->fb_type);
    MPI_Comm_free (&configMPI->cart);
}

void calcInnerGrid(heatGrid* grid, dataMPI* mympi, double dt)
{
    double dtheta;

    /* Calculate time step for inner cells: read from theta, write new timestep to thetanew. Dont calc outer cells. */
    for(int x = 2; x < grid->local_size_x; x++) {
        for(int y = 2; y < grid->local_size_y; y++) {
          for(int z = 2; z < grid->local_size_z; z++) {
            dtheta = ( grid->theta[x-1][y][z] + grid->theta[x+1][y][z] - 2*grid->theta[x][y][z]) / (grid->dx * grid->dx)
                   + ( grid->theta[x][y-1][z] + grid->theta[x][y+1][z] - 2*grid->theta[x][y][z]) / (grid->dy * grid->dy)
                   + ( grid->theta[x][y][z-1] + grid->theta[x][y][z+1] - 2*grid->theta[x][y][z]) / (grid->dz * grid->dz);
            grid->thetanew[x][y][z] = grid->theta[x][y][z] + grid->k * dtheta * dt;
        }
      }
    }
}

// calc inner cells of the row at position row, dont calc corners because they depend on other ghost cells
void calc_tb(heatGrid* grid, dataMPI* mympi, double dt, int tb)
{
    double dtheta;

    for(int y = 1; y < grid->local_size_y + 1; ++y) {
      for(int z = 1; z < grid->local_size_z + 1; ++z) {
        dtheta = (grid->theta[tb-1][y][z] + grid->theta[tb+1][y][z] - 2*grid->theta[tb][y][z]) / (grid->dx * grid->dx)
               + (grid->theta[tb][y-1][z] + grid->theta[tb][y+1][z] - 2*grid->theta[tb][y][z]) / (grid->dy * grid->dy)
               + (grid->theta[tb][y][z-1] + grid->theta[tb][y][z+1] - 2*grid->theta[tb][y][z]) / (grid->dz * grid->dz);
        grid->thetanew[tb][y][z] = grid->theta[tb][y][z] + grid->k * dtheta * dt;
      }
    }
}

void calc_lr(heatGrid* grid, dataMPI* mympi, double dt, int lr)
{
    double dtheta;

    for(int x = 2; x < grid->local_size_x; ++x) {
      for(int z = 1; z < grid->local_size_z + 1; ++z) {
        dtheta = (grid->theta[x-1][lr][z] + grid->theta[x+1][lr][z] - 2*grid->theta[x][lr][z]) / (grid->dx * grid->dx)
               + (grid->theta[x][lr-1][z] + grid->theta[x][lr+1][z] - 2*grid->theta[x][lr][z]) / (grid->dy * grid->dy)
               + (grid->theta[x][lr][z-1] + grid->theta[x][lr][z+1] - 2*grid->theta[x][lr][z]) / (grid->dz * grid->dz);
        grid->thetanew[x][lr][z] = grid->theta[x][lr][z] + grid->k * dtheta * dt;
      }
    }
}

void calc_fb(heatGrid* grid, dataMPI* mympi, double dt, int fb)
{

    double dtheta;

    for(int x = 2; x < grid->local_size_x; ++x) {
      for(int y = 2; y < grid->local_size_y; ++y) {
        dtheta = (grid->theta[x-1][y][fb] + grid->theta[x+1][y][fb] - 2*grid->theta[x][y][fb]) / (grid->dx * grid->dx)
               + (grid->theta[x][y-1][fb] + grid->theta[x][y+1][fb] - 2*grid->theta[x][y][fb]) / (grid->dy * grid->dy)
               + (grid->theta[x][y][fb-1] + grid->theta[x][y][fb+1] - 2*grid->theta[x][y][fb]) / (grid->dz * grid->dz);
        grid->thetanew[x][y][fb] = grid->theta[x][y][fb] + grid->k * dtheta * dt;
      }
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

    if(argc != 3) {
        perror("Usage: heatMPI <size_per_dim> <iterations>");
        return 1;
    }

    mygrid.global_size_x = strtol(argv[1], NULL, 10);
    mygrid.global_size_y = strtol(argv[1], NULL, 10);
    mygrid.global_size_z = strtol(argv[1], NULL, 10);

    /* setup MPI */
    heatMPISetup(&mygrid, &mympi);

    /* create heatGrid and initialize variables */
    heatAllocate(&mygrid);
    nsteps = strtol(argv[2], NULL, 10);
    energyInitial = energyFinal = 0.0;
    dt = 0.05;

    if (mympi.rank == 0)
        heatInitialize(&mygrid);

    /* energy of initial grid */
    heatTotalEnergy(&mygrid, &energyLocal);
    MPI_Reduce(&energyLocal, &energyInitial, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (mympi.rank == 0) {

        printf("Running Jacobi on %d processors with (%d, %d %d) elements.\n", mympi.size, mympi.num_cells_x, mympi.num_cells_y, mympi.num_cells_z);
        printf("Array Dimensions: %d %d %d\n", mygrid.global_size_x, mygrid.global_size_y, mygrid.global_size_z);
        printf("Block Dimensions: %d %d %d\n", mygrid.local_size_x, mygrid.local_size_y, mygrid.local_size_z);
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
        MPI_Request send_requests[6];
        MPI_Request recv_requests[6];

        //Irecv before Isend for better performance

        /* top and bottom */

        /*Receive Left border column from left neighbor*/
        MPI_Irecv(&(mygrid.theta[0][1][1]),
                1, mympi.tb_type, mympi.up, step, mympi.cart, &recv_requests[0]);
        /*Send right column to right neighbor*/
        MPI_Isend(&(mygrid.theta[mygrid.local_size_x][1][1]),
                1, mympi.tb_type, mympi.down, step, mympi.cart, &send_requests[0]);

        /*Receive Left border column from left neighbor*/
        MPI_Irecv(&(mygrid.theta[mygrid.local_size_x+1][1][1]),
                1, mympi.tb_type, mympi.down, step, mympi.cart, &recv_requests[1]);
        /*Send right column to right neighbor*/
        MPI_Isend(&(mygrid.theta[1][1][1]),
                1, mympi.tb_type, mympi.up, step, mympi.cart, &send_requests[1]);

        /* left and right */
        MPI_Irecv(&(mygrid.theta[1][0][1]),
                1, mympi.lr_type, mympi.left, step, mympi.cart, &recv_requests[2]);
        /*Send right column to right neighbor*/
        MPI_Isend(&(mygrid.theta[1][mygrid.local_size_y][1]),
                1, mympi.lr_type,  mympi.right, step, mympi.cart, &send_requests[2]);

        MPI_Irecv(&(mygrid.theta[1][mygrid.local_size_y+1][1]),
                1, mympi.lr_type, mympi.right, step, mympi.cart, &recv_requests[3]);
        /*Send right column to right neighbor*/
        MPI_Isend(&(mygrid.theta[1][1][1]),
                1, mympi.lr_type,  mympi.left, step, mympi.cart, &send_requests[3]);

        /* front and back */

        MPI_Irecv(&(mygrid.theta[1][1][0]),
                1, mympi.fb_type, mympi.back, step, mympi.cart, &recv_requests[4]);
        /*Send right column to right neighbor*/
        MPI_Isend(&(mygrid.theta[1][1][mygrid.local_size_z]),
                1, mympi.fb_type,  mympi.front, step, mympi.cart, &send_requests[4]);

        MPI_Irecv(&(mygrid.theta[1][1][mygrid.local_size_z+1]),
                1, mympi.fb_type, mympi.front, step, mympi.cart, &recv_requests[5]);
        /*Send right column to right neighbor*/
        MPI_Isend(&(mygrid.theta[1][1][1]),
                1, mympi.fb_type,  mympi.back, step, mympi.cart, &send_requests[5]);

        minimon.leave("calc update async");
        //overlap communication with computation
        minimon.enter();
        calcInnerGrid(&mygrid, &mympi, dt);
        minimon.leave("calc inner");

        minimon.enter();
        // wait for all Recv to complete and save switch-case statement
        MPI_Waitall(6, recv_requests, MPI_STATUSES_IGNORE);
        minimon.leave("calc wait");

        /* calcColumn has more performance than calcRow because of memory access pattern,
           thus using calcColumn to calc corners */

        minimon.enter();
        calc_tb(&mygrid, &mympi, dt,  1);
        calc_tb(&mygrid, &mympi, dt,  mygrid.local_size_x);
        calc_lr(&mygrid, &mympi, dt,  1);
        calc_lr(&mygrid, &mympi, dt,  mygrid.local_size_y);
        calc_fb(&mygrid, &mympi, dt,  1);
        calc_fb(&mygrid, &mympi, dt,  mygrid.local_size_z);
        minimon.leave("calc bound");

        //Sends have to be finished before swaping pointers
        MPI_Waitall(6, send_requests, MPI_STATUSES_IGNORE);

        //swap pointer
        double ***tmp;
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
