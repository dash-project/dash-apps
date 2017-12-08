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

set output "trace.png"
#set title "Foo"

set style data points
set datafile separator ";"
set logscale x 10
set logscale y 10
set key right bottom Left

f(x,y)= (x eq "P1"? y : 1/0)

plot \
    'trace.csv' using 7:( strcol(2) eq "main" ? $6 : 1/0 ) with points ps 1 t "main", \
    'trace.csv' using 7:( strcol(2) eq "dash::init" ? $6 : 1/0 ) with points ps 1 t "dash::init", \
    'trace.csv' using 7:( strcol(2) eq "dash::finalize" ? $6 : 1/0 ) with points ps 1 t "dash::final", \
    'trace.csv' using 7:( strcol(2) eq "scaledown" ? $6 : 1/0 ) with points ps 1 t "scaledown", \
    'trace.csv' using 7:( strcol(2) eq "scaleup" ? $6 : 1/0 ) with points ps 1 t "scaleup", \
    'trace.csv' using 7:( strcol(2) eq "setup" ? $6 : 1/0 ) with points ps 1 t "setup", \
    'trace.csv' using 7:( strcol(2) eq "smooth_col_bc" ? $6 : 1/0 ) with points ps 1 t "smooth col bc", \
    'trace.csv' using 7:( strcol(2) eq "smooth_inner" ? $6 : 1/0 ) with points ps 1 t "smooth inne", \
    'trace.csv' using 7:( strcol(2) eq "smooth_outer" ? $6 : 1/0 ) with points ps 1 t "smooth oute", \
    'trace.csv' using 7:( strcol(2) eq "smooth_wait" ? $6 : 1/0 ) with points ps 1 t "smooth wait", \
    'trace.csv' using 7:( strcol(2) eq "smoothen" ? $6 : 1/0 ) with points ps 1 t "smoothen", \
    'trace.csv' using 7:( strcol(2) eq "smoothfinal" ? $6 : 1/0 ) with points ps 1 t "smooth final"

set output "trace_sep.png"
plot \
    'trace.csv' using 7:( strcol(2) eq "smooth_inner" ? $6 : 1/0 ) with points ps 1 t "smooth inne", \
    'trace.csv' using 7:( strcol(2) eq "smooth_outer" ? $6 : 1/0 ) with points ps 1 t "smooth oute", \
    'trace.csv' using 7:( strcol(2) eq "smooth_wait" ? $6 : 1/0 ) with points ps 1 t "smooth wait", \
    'trace.csv' using 7:( strcol(2) eq "smoothen" ? $6 : 1/0 ) with points ps 1 t "smoothen", \

set ylabel "Flop/s"

set output "trace_flops.png"
plot \
    'trace.csv' using 7:( strcol(2) eq "smooth_inner" ? $8/$6 : 1/0 ) with points ps 1 t "smooth inne", \
    'trace.csv' using 7:( strcol(2) eq "smooth_outer" ? $8/$6 : 1/0 ) with points ps 1 t "smooth oute", \
    'trace.csv' using 7:( strcol(2) eq "smoothen" ? $8/$6 : 1/0 ) with points ps 1 t "smoothen", \
