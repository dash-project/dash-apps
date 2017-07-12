To view the 3D cube download paraview from https://www.paraview.org/download/

After executing multigrid3d open them from paraview with File -> open. Navigate to your build folder
and select the image group as a whole. 

In the subwindow "Pipeline Browser" select image.csv* and then press Apply ins the Properties window below.

Now add a filter to the selected group with: Filters -> Alphabetical -> Table to Points.

In the "Properties" window specify the x,y and z collumns with the corresponding collumns int the csv and press Apply.

If the cube does not show up in the Layout #1 screen, select it and click the "eye"-icon next to the TableToPoints filter.

In the "Properties" window scroll down and in the section "Coloring" change Solid Color to heat.
