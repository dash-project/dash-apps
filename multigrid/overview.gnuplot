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
    filename using 4:( strcol(2) eq "smoothen" ? $7 : 1/0 ):8:9 with yerrorbars t "smoothen", \
    filename using 4:( strcol(2) eq "smoothen_inner" ? $7 : 1/0 ):8:9 with yerrorbars t "smoothen inner", \
    filename using 4:( strcol(2) eq "smoothen_outer" ? $7 : 1/0 ):8:9 with yerrorbars t "smoothen outer", \
    filename using 4:( strcol(2) eq "scaledown" ? $7 : 1/0 ):8:9 with yerrorbars t "scaledown", \
    filename using 4:( strcol(2) eq "scaleup" ? $7 : 1/0 ):8:9 with yerrorbars t "scaleup", \
    filename using 4:( strcol(2) eq "main" ? $7 : 1/0 ):8:9 with yerrorbars t "main", \
    filename using 4:( strcol(2) eq "setup" ? $7 : 1/0 ):8:9 with yerrorbars t "setup", \
    filename using 4:( strcol(2) eq "dash::init" ? $7 : 1/0 ):8:9 with yerrorbars t "dash::init", \
    filename using 4:( strcol(2) eq "dash::finalize" ? $7 : 1/0 ):8:9 with yerrorbars t "dash::finalize"


filename_out=filename.".detail.png"
set output filename_out
#set title "Foo"

set xrange[0.5:*]
plot \
    filename using 4:( strcol(2) eq "smoothen" ? $7 : 1/0 ):8:9 with yerrorbars t "smoothen", \
    filename using 4:( strcol(2) eq "smoothen_inner" ? $7 : 1/0 ):8:9 with yerrorbars t "smoothen inner", \
    filename using 4:( strcol(2) eq "smoothen_outer" ? $7 : 1/0 ):8:9 with yerrorbars t "smoothen outer", \
    filename using 4:( strcol(2) eq "smoothen_wait" ? $7 : 1/0 ):8:9 with yerrorbars t "smoothen wait", \
    filename using 4:( strcol(2) eq "smoothen_wait_res" ? $7 : 1/0 ):8:9 with yerrorbars t "smoothen wait res", \
    filename using 4:( strcol(2) eq "smoothen_collect" ? $7 : 1/0 ):8:9 with yerrorbars t "smoothen collect", \
    filename using 4:( strcol(2) eq "scaledown" ? $7 : 1/0 ):8:9 with yerrorbars t "scaledown", \
    filename using 4:( strcol(2) eq "scaleup" ? $7 : 1/0 ):8:9 with yerrorbars t "scaleup"

filename_out=filename.".flops.png"
set output filename_out
#set title "Foo"
set ylabel "Speed [Flop/s]"

set xrange[0.5:*]
plot \
    filename using 3:( strcol(1) eq "smoothen" ? $5/$7 : 1/0 ):($5/$8):($5/$9) with yerrorbars t "smoothen", \
    filename using 3:( strcol(1) eq "smoothen_inner" ? $5/$7 : 1/0 ):($5/$8):($5/$9) with yerrorbars t "smoothen inner", \
    filename using 3:( strcol(1) eq "smoothen_outer" ? $5/$7 : 1/0 ):($5/$8):($5/$9) with yerrorbars t "smoothen outer", \
    filename using 3:( strcol(1) eq "scaledown" ? $5/$7 : 1/0 ):($5/$8):($5/$9) with yerrorbars t "scaledown", \
    filename using 3:( strcol(1) eq "scaleup" ? $5/$7 : 1/0 ):($5/$8):($5/$9) with yerrorbars t "scaleup"












