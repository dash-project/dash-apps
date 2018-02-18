#!/usr/bin/gnuplot

# To give the input file name on the command line call it like this:
# gnuplot -e "filename='<filename>'" overview.gnuplot


set terminal png size 1200,900


set xlabel "# grid elements"
set ylabel "time [s]"
set grid

if (!exists("filename")) filename='overview.csv'


set datafile separator ";"
set logscale xy 10
set key left top Left


filename_out=filename.".png"
set output filename_out
#set title "Foo"

set xrange[0.5:*]
plot \
    filename using 3:( strcol(1) eq "smoothen" ? $6 : 1/0 ):7:8 with yerrorbars t "smoothen", \
    filename using 3:( strcol(1) eq "smoothen_inner" ? $6 : 1/0 ):7:8 with yerrorbars t "smoothen inner", \
    filename using 3:( strcol(1) eq "smoothen_outer" ? $6 : 1/0 ):7:8 with yerrorbars t "smoothen outer", \
    filename using 3:( strcol(1) eq "scaledown" ? $6 : 1/0 ):7:8 with yerrorbars t "scaledown", \
    filename using 3:( strcol(1) eq "scaleup" ? $6 : 1/0 ):7:8 with yerrorbars t "scaleup", \
    filename using 3:( strcol(1) eq "main" ? $6 : 1/0 ):7:8 with yerrorbars t "main", \
    filename using 3:( strcol(1) eq "setup" ? $6 : 1/0 ):7:8 with yerrorbars t "setup", \
    filename using 3:( strcol(1) eq "dash::init" ? $6 : 1/0 ):7:8 with yerrorbars t "dash::init", \
    filename using 3:( strcol(1) eq "dash::finalize" ? $6 : 1/0 ):7:8 with yerrorbars t "dash::finalize"


filename_out=filename.".detail.png"
set output filename_out
#set title "Foo"

set xrange[0.5:*]
plot \
    filename using 3:( strcol(1) eq "smoothen" ? $6 : 1/0 ):7:8 with yerrorbars t "smoothen", \
    filename using 3:( strcol(1) eq "smoothen_inner" ? $6 : 1/0 ):7:8 with yerrorbars t "smoothen inner", \
    filename using 3:( strcol(1) eq "smoothen_outer" ? $6 : 1/0 ):7:8 with yerrorbars t "smoothen outer", \
    filename using 3:( strcol(1) eq "smoothen_wait" ? $6 : 1/0 ):7:8 with yerrorbars t "smoothen wait", \
    filename using 3:( strcol(1) eq "smoothen_wait_res" ? $6 : 1/0 ):7:8 with yerrorbars t "smoothen wait res", \
    filename using 3:( strcol(1) eq "smoothen_collect" ? $6 : 1/0 ):7:8 with yerrorbars t "smoothen collect", \
    filename using 3:( strcol(1) eq "scaledown" ? $6 : 1/0 ):7:8 with yerrorbars t "scaledown", \
    filename using 3:( strcol(1) eq "scaleup" ? $6 : 1/0 ):7:8 with yerrorbars t "scaleup"

filename_out=filename.".flops.png"
set output filename_out
#set title "Foo"

set xrange[0.5:*]
plot \
    filename using 3:( strcol(1) eq "smoothen" ? $4/$6 : 1/0 ):($4/$7):($4/$8) with yerrorbars t "smoothen", \
    filename using 3:( strcol(1) eq "smoothen_inner" ? $4/$6 : 1/0 ):($4/$7):($4/$8) with yerrorbars t "smoothen inner", \
    filename using 3:( strcol(1) eq "smoothen_outer" ? $4/$6 : 1/0 ):($4/$7):($4/$8) with yerrorbars t "smoothen outer", \
    filename using 3:( strcol(1) eq "scaledown" ? $4/$6 : 1/0 ):($4/$7):($4/$8) with yerrorbars t "scaledown", \
    filename using 3:( strcol(1) eq "scaleup" ? $4/$6 : 1/0 ):($4/$7):($4/$8) with yerrorbars t "scaleup"












