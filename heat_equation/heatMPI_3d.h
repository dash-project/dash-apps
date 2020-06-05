#include <math.h>
#include <string.h>
#include <stdio.h>
#include <malloc.h>
#include <mpi.h>

/* This type represents the grid and its description */
typedef struct {
    /* current theta array */
    double ***theta;
    /* new theta array */
    double ***thetanew;
    /* total size of all arrays of all processes combined */
    int global_size_x;
    int global_size_y;
    int global_size_z;
    /* actual local size ot the array managed by one process without ghost cells */
    int local_size_x;
    int local_size_y;
    int local_size_z;
    /* size of a grid cell */
    double dx;
    double dy;
    double dz;
    /* "heat equation constant" */
    double k;
} heatGrid;

/* This type defines everything thats needed for MPI */
typedef struct {
    /* Own rank, used to only let master do output*/
    int rank;
    /* Number of all processes*/
    int size;
    /* Comm for a cartesian distribution of the grid*/
    MPI_Comm cart;
    /* Neighbors in communicator*/
    int up,down,left,right, front, back;
    /* Dimensions of cartesian array*/
    int num_cells_x;
    int num_cells_y;
    int num_cells_z;

    int dim_x;
    int dim_y;
    int dim_z;
    /* Position of process in cartesian array (x, y)*/
    int coord_x;
    int coord_y;
    int coord_z;
    /* Datatype used to transfer a data column*/
    MPI_Datatype tb_type; // top bottom (tb)
    MPI_Datatype lr_type; // left right (lr)
    MPI_Datatype fb_type; // front back (fb)
} dataMPI;

void heatAllocate(heatGrid* grid);
void heatDeallocate(heatGrid* grid);

void heatInitialize(heatGrid* grid);

void heatPrint(heatGrid* grid);
void heatTotalEnergy(heatGrid* grid, double *energy);

void heatMPISetup (heatGrid* grid, dataMPI* configMPI);
void heatMPIFree (dataMPI* configMPI);

void calcInnerGrid(heatGrid* grid, dataMPI* mympi, double dt);

void calc_tb(heatGrid* grid, dataMPI* mympi, double dt, int tb);
void calc_lr(heatGrid* grid, dataMPI* mympi, double dt, int lr);
void calc_fb(heatGrid* grid, dataMPI* mympi, double dt, int fb);
