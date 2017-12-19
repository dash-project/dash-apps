#!/usr/bin/gnuplot

set terminal png size 1200,900


set xlabel "# grid elements"
set ylabel "time [s]"
set grid

#set style line idx_c     linecolor rgb "#ff0000" linetype 0
#set style line idx_c_ikj linecolor rgb "#ff0000" linetype 1
#set style line idx_f     linecolor rgb "#00ff00" linetype 0
#set style line idx_f_jki linecolor rgb "#00ff00" linetype 1
#set style line idx_j     linecolor rgb "#0000ff" linetype 1
#set style line idx_h     linecolor rgb "#00ffff" linetype 1
if (!exists("filename")) filename='trace.csv'

filename_out=filename.".png"
set output filename_out
#set title "Foo"

set style data points
set datafile separator ";"
set logscale x 10
set logscale y 10
set key right bottom Left

#f(x,y)= (x eq "P1"? y : 1/0)

plot filename using 5:( strcol(2) eq "main" ? $4 : 1/0 ) with points ps 1 t "main", \
    filename using 5:( strcol(2) eq "dash::init" ? $4 : 1/0 ) with points ps 1 t "dash::init", \
    filename using 5:( strcol(2) eq "dash::final" ? $4 : 1/0 ) with points ps 1 t "dash::final", \
    filename using 5:( strcol(2) eq "scaledown" ? $4 : 1/0 ) with points ps 1 t "scaledown", \
    filename using 5:( strcol(2) eq "scaleup" ? $4 : 1/0 ) with points ps 1 t "scaleup", \
    filename using 5:( strcol(2) eq "setup" ? $4 : 1/0 ) with points ps 1 t "setup", \
    filename using 5:( strcol(2) eq "smooth_col_bc" ? $4 : 1/0 ) with points ps 1 t "smooth col bc", \
    filename using 5:( strcol(2) eq "smooth_inner" ? $4 : 1/0 ) with points ps 1 t "smooth inner", \
    filename using 5:( strcol(2) eq "smooth_outer" ? $4 : 1/0 ) with points ps 1 t "smooth outer", \
    filename using 5:( strcol(2) eq "smooth_wait" ? $4 : 1/0 ) with points ps 1 t "smooth wait", \
    filename using 5:( strcol(2) eq "smoothen" ? $4 : 1/0 ) with points ps 1 t "smoothen", \
    filename using 5:( strcol(2) eq "smoothen_res" ? $4 : 1/0 ) with points ps 1 t "smoothen res", \
    filename using 5:( strcol(2) eq "smoothfinal" ? $4 : 1/0 ) with points ps 1 t "smooth final"

filename_out=filename."_seq.png"
set output filename_out
plot \
    filename using 5:( strcol(2) eq "smooth_inner" ? $4 : 1/0 ) with points ps 1 t "smooth inner", \
    filename using 5:( strcol(2) eq "smooth_outer" ? $4 : 1/0 ) with points ps 1 t "smooth outer", \
    filename using 5:( strcol(2) eq "smooth_wait" ? $4 : 1/0 ) with points ps 1 t "smooth wait", \
    filename using 5:( strcol(2) eq "smoothen" ? $4 : 1/0 ) with points ps 1 t "smoothen", \
    filename using 5:( strcol(2) eq "smoothen res" ? $4 : 1/0 ) with points ps 1 t "smoothen res"

set ylabel "Flop/s"

filename_out=filename."_flops.png"
set output filename_out
plot \
    filename using 5:( strcol(2) eq "smooth_inner" ? $6/$4 : 1/0 ) with points ps 1 t "smooth inner", \
    filename using 5:( strcol(2) eq "smooth_outer" ? $6/$4 : 1/0 ) with points ps 1 t "smooth outer", \
    filename using 5:( strcol(2) eq "smoothen" ? $6/$4 : 1/0 ) with points ps 1 t "smoothen", \
    filename using 5:( strcol(2) eq "smoothen res" ? $6/$4 : 1/0 ) with points ps 1 t "smoothen res"

filename_out=filename."_inner.png"
set output filename_out
plot filename using 5:( strcol(2) eq "smooth_inner" ? $4 : 1/0 ) with points ps 1 t "smooth inner"

filename_out=filename."_smoothen.png"
set output filename_out
plot filename using 5:( strcol(2) eq "smoothen" ? $4 : 1/0 ) with points ps 1 t "smoothen"


