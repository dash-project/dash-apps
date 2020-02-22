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
 *  Header for MPI version.
 */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <malloc.h>
#include <mpi.h>

/* This type represents the grid and its description */
typedef struct {
    /* current theta array */
    double **theta;
    /* new theta array */
    double **thetanew;
    /* total size of all arrays of all processes combined */
    int global_size_x;
    int global_size_y;
    /* actual local size ot the array managed by one process without ghost cells */
    int local_size_x;
    int local_size_y;
    /* size of a grid cell */
    double dx;
    double dy;
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
    int up,down,left,right;
    /* Dimensions of cartesian array*/
    int num_cells_x;
    int num_cells_y;

    int dim_x;
    int dim_y;
    /* Position of process in cartesian array (x, y)*/
    int coord_x;
    int coord_y;
    /* Datatype used to transfer a data column*/
    MPI_Datatype rowtype;
} dataMPI;

void heatAllocate(heatGrid* grid);
void heatDeallocate(heatGrid* grid);

void heatInitialize(heatGrid* grid);

void heatPrint(heatGrid* grid);
void heatTotalEnergy(heatGrid* grid, double *energy);

void heatMPISetup (heatGrid* grid, dataMPI* configMPI);
void heatMPIFree (dataMPI* configMPI);

void calcInnerGrid(heatGrid* grid, dataMPI* mympi, double dt);
void calcRow(heatGrid* grid, dataMPI* mympi, double dt, int row);
void calcColumn(heatGrid* grid, dataMPI* mympi, double dt, int column, int start, int end);
void calcCorners(heatGrid* grid, dataMPI* mympi, double dt);
