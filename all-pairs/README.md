# All-Pairs
Measures the latency between each pair of MPI ranks using various kernels.
The results are stored as hdf5 files.
To analyze the data, use the provided r-script `load_and_process.r`

## Build
This project uses cmake to setup the build environment. It depends on Boost (program_options, log) and a build of DASH-MPI which can be downloaded here: http://dash-project.org/
When everything is set up correctly, run build.sh to build the project using cmake.

### Troubleshooting
If your version of Boost is too new, cmake might not detect the dependencies correctly. Hence the application does not link correctly.
In this case, uncomment the provided workaround in `CMakeLists.txt`.

## Run
For the collection of the data run the application with `mpiexec ./all-pairs <params>`. The output files are placed in the current directory.
As this application uses parallel I/O, run it on a parallel filesystem only.

For the postprocessing see the README in the r-scripts folder.

### Parameters
To show all options, run `all-pairs --help`

## Add Kernels
Feel free to write your own kernels by extending `AllPairsKernel` or a subclass of it.

## Examples
### LRZ SuperMUC Phase II HPC system
Run on 20 haswell nodes of [SuperMUC](https://www.lrz.de/services/compute/supermuc/systemdescription/), using intel MPI. The cassis switches (3) can be easily seen as light blue blocks.

![supermuc 20 nodes](https://gist.githubusercontent.com/fmoessbauer/051bc8340e0d6c05432394b95f521bc2/raw/414acc1f12c0acf7551224230f5746a65d7f4f46/supermuc_20_nodes_RMA_GET_web.png)

### HLRS Hazel Hen HPC system
Run on 5 nodes of [Hazel Hen](https://www.hlrs.de/de/systems/cray-xc40-hazel-hen/) using cray MPI. The schedule is disadvantageous because the nodes are split over two groups of the dragonfly interconnect.

![hazelhen 5 nodes]
(https://gist.githubusercontent.com/fmoessbauer/051bc8340e0d6c05432394b95f521bc2/raw/414acc1f12c0acf7551224230f5746a65d7f4f46/hazelhen_5_nodes_RMA_GET_web.png)
