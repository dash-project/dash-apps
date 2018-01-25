# How to view 3D views of the multigrid 3D example

* Install paraview from distro or download paraview from https://www.paraview.org/download/

* After executing multigrid3d open them from paraview with File -> open. Navigate to your build folder and select the image group as a whole aas 'image.csv'.

* In the subwindow "Pipeline Browser" select image.csv* and then press Apply ins the Properties window below.

* Now add a filter to the selected group with: Filters -> Alphabetical -> Table to Structured Grid.

* In the "Properties" window specify the x,y and z columns with the corresponding columns in the csv. Also specify the correct extents in x-, y-, and z- direction from 0 to 2^n-1 -- you'll get an error message if it is incorrect. Then switch Representation to "Outline" and press Apply.

* In the "Properties" window scroll down and in the section "Coloring" change Solid Color to heat.

* If the cube does not show up in the Layout #1 screen, select it and click the "eye"-icon next to the entries in the Pipeline Browser.

* Now add a filter to the selected group with: Filters -> Alphabetical -> Iso Volume. Then Set Minimum to 6 and Maximum to 10 -- that will be the two isosurfaces that you'll see. Farther down set Opacity to 0.5.

* Afterwards play with the green time forward/backward arrows to step through the evolution of the grids on every level. You can also "Save Animation".

# How to visualize the minimon trace data with gnuplot

* Use the script combine_csvs.sh to merge all trace files into one
  usage: ./combine_csvs.sh <name extension>

* Build all plots with gnuplot
  usage: gnuplot -e "filename='<name_of_tracefile.csv>'" trace.gnuplot
