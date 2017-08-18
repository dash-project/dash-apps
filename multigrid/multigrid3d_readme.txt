# How to view 3D views of the multigrid 3D example

* Install paraview from distro or download paraview from https://www.paraview.org/download/

* After executing multigrid3d open them from paraview with File -> open. Navigate to your build folder and select the image group as a whole aas 'image.csv'.

* In the subwindow "Pipeline Browser" select image.csv* and then press Apply ins the Properties window below.

* Now add a filter to the selected group with: Filters -> Alphabetical -> Table to Points.

* In the "Properties" window specify the x,y and z columns with the corresponding columns in the csv. Then switch Representation to "Outline". Then and press Apply.

* In the "Properties" window scroll down and in the section "Coloring" change Solid Color to heat.

* If the cube does not show up in the Layout #1 screen, select it and click the "eye"-icon next to the entries in the Pipeline Browser.

* Now add a filter to the selected group with: Filters -> Alphabetical -> Iso Volume. Then Set Minimum to 6 and Maximum to 10 and farther down Opacity to 0.5.

* Afterwards play with the green time forward/backward arrows to step through the evolution of the grids on every level.
