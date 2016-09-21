# All-Pairs
Measures the latency between each pair of MPI ranks using various kernels.
The results are stored as hdf5 files.
To analyze the data, use the provided r-script `load_and_process.r`

## Build
This project uses cmake to setup the build environment. It depends on a build of DASH-MPI which can be downloaded here: http://dash-project.org/
When everything is set up correctly, run build.sh to build the project using cmake.

## Run
For the collection of the data run the application with mpiexec `./all-pairs <params>`. The output files are placed in the current directory.
As this application uses parallel I/O, run it on a parallel filesystem only.

For the postprocessing see the README in the r-scripts folder.

### Parameters
To show all options, run `all-pairs --help`

## Add Kernels
Feel free to write your own kernels by extending `AllPairsKernel` or a subclass of it.
