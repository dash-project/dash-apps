# The multigrid 3d example code

This DASH example code solves the heat equation in a 3D cubiod. It either uses
* the multigrid method to compute the steady state
* a iterative solver to compute the steady state
* a simulation mode doing a time simulation

## Command line options:

     Call me as 'mpirun ./multigrid3d' [-h|--help] [levels(default 5)] [...more options...]
     <l>           number of levels to use at most, the simulation grid which is
                   the finest grid in the multigrid hierarchy of grids will have
                   2^l -1 inner elements plus 2 boundary elements per dimension

     Modes of operation

     -t|--test     run some internal tests
     -e|--elastic  use elastic multigrid mode i.e., use fewer units (processes)
                   on coarser grids
     -f|--flat     run flat mode i.e., use iterative solver on a single grid
     --sim <t> <s> run a simulation over time, that is also a "flat" solver
                   working only on a single grid. It runs t seconds simulation
                   time. The time step dt is determined by the grid and the
                   stability condition. This mode matches all time steps n*s <= t
                   exactly for the sake of a nice visualization.
                   (Visualization only active when compiled with WITHCSVOUTPUT.)

     Further options

     --eps <eps>   define epsilon for the iterative solver in flat or multigrid modes,
                   the iterative solver on any grid stops when residual <= eps
     -g <n>        determine size of CSV output -- only when compiled with WITHCSVOUTPUT
                   the combined CVS output from all units (processes) will be
                   2^n +1 points in every dimension, default is n= 5 or 33^3 elements
     -d <d h w>    Set physical dimensions of the simulation grid in meters
                   (default 10.0, 10.0, 10.0)

     (This executable was compiled without WITHCSVOUTPUT)

## Simulation mode with '--sim'

In simulation mode the heat equation is simulated for t seconds with an output every s seconds (if compiled with WITHCSVOUTPUT). The time step in seconds is adjusted to the square of the spacial distance h automatically according to the stability criterion.

The simulation mode works on a single grid, not multigrid, not using a hierarchy of grids.

### I/O behavior in simulation mode

I simulation mode with WITHCSVOUTPUT enabled there is a regular I/O pattern where every unit writes a separate file per output time step, which can be used for visualization with Paraview -- see below.

* The 'l' option determines the global grid size as l --> (2^l -1)³ inner grid elements (in 3 dimensions) or (2^l +1)³ total grid elements including the boundary values (one each at the lower and upper end per dimension). This defines the effort for every smoothing step on the grid. (Hint: Increase l to have notably larger gaps between the I/O phases due to much larger computation. It also requires much more memory, though.)
The application can use 2^(l-2) units (processes) at max. However, it gets rather inefficient when getting close to it. Recommended is something like 2^(l-9) to 2^(l-7) units, that is 513³ to 129³ grid elements per unit which requires approx. 3 GB to 48 MB per unit.

* The '-d' option defines the physical dimensions of the simulation area. This is not the number of grid elements per dimension as defined by the 'l' argument (see previous item). But the physical dimension influences the grid distance h and h² determines the time step dt (delta t). Thus smaller physical dimensions for the same compute grid cause much smaller time steps. Thus, more solver iterations are required to cover the same simulation time period until the next output time step. If one (actually the shortest) of the physical dimensions is divided by two then the time step dt is divided by 4. (Hint: Keep the default.)

* The '-g' option determines the global output grid size. In the code the output grid size can be arbitrary. With the '-g' option it can be set to (2^n+1)³ by specifying 'n'.   Default is n= 5 or 33³ elements. At every output time step every of the p units (processes) write (2^n+1)³/p elements. It write one ASCII line per element with 24 to 42 bytes in a single file write operation. (Hint: Increase the argument 'n' to produce much larger I/O operations per process. That is 8-fold when increasing n by 1.)

* Finally, the '--sim' option has two arguments t and s. With 't' one determines the total simulation time period in simulated seconds. With s one specifies the time between two output time steps e.g., with s=0.04 one gets 25 FPS for a video with real time. ... producing rather boring videos ;)
(Hint: Increase t for longer runs, use s to closely control the period between to I/O phases.)

## How to produce output for visualization

When compiled with '-DWITHCSVOUTPUT' (or with the given Makefile using the multigrid3d_csv executable) the program writes out CSV output for visualization.
* Every unit (process) writes out the local part of the simulation grid including the boundary values to a separate file.
* The 'combine_csvs.sh' script combines the per-unit local parts for a time step into a single file per time step.
* Paraview can read this as a single file or as many successive time steps -- see below.


## How to view 3D views of the multigrid 3D example

* Install paraview from distro or download paraview from https://www.paraview.org/download/

* After executing multigrid3d open them from paraview with File -> open. Navigate to your build folder and select the image group as a whole as 'image.csv'.

* In the subwindow "Pipeline Browser" select image.csv* and then press Apply in the Properties window below.

* Now add a filter to the selected group with: Filters -> Alphabetical -> Table to Structured Grid.

* In the "Properties" window specify the x,y and z columns with the corresponding columns in the csv. Also specify the correct extents in x-, y-, and z- direction from 0 to 2^n-1 -- you'll get an error message if it is incorrect. Then switch Representation to "Outline" and press Apply.

* In the "Properties" window scroll down and in the section "Coloring" change Solid Color to heat.

* If the cube does not show up in the Layout #1 screen, select it and click the "eye"-icon next to the entries in the Pipeline Browser.

* Now add a filter to the selected group with: Filters -> Alphabetical -> Iso Volume. Then Set Minimum to 3 and Maximum to 7 -- that will be the two isosurfaces that you'll see. Farther down set Opacity to 0.5.

* Afterwards play with the green time forward/backward arrows to step through the evolution of the grids on every level. You can also "Save Animation".

## How to visualize the minimon trace data with gnuplot

Use the script combine_csvs.sh (same as above) to merge all trace files into one and have gnuplot create PNG output files

Usage:

    ./combine_csvs.sh <name extension>
