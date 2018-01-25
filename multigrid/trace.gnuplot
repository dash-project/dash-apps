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
set logscale xy 10
set key left top Left

#f(x,y)= (x eq "P1"? y : 1/0)

plot \
    filename using 5:( strcol(2) eq "main" ? $4 : 1/0 ) with points ps 1 t "main", \
    filename using 5:( strcol(2) eq "dash::init" ? $4 : 1/0 ) with points ps 1 t "dash::init", \
    filename using 5:( strcol(2) eq "dash::final" ? $4 : 1/0 ) with points ps 1 t "dash::final", \
    filename using 5:( strcol(2) eq "scaledown" ? $4 : 1/0 ) with points ps 1 t "scaledown", \
    filename using 5:( strcol(2) eq "scaleup" ? $4 : 1/0 ) with points ps 1 t "scaleup", \
    filename using 5:( strcol(2) eq "setup" ? $4 : 1/0 ) with points ps 1 t "setup", \
    filename using 5:( strcol(2) eq "smooth" ? $4 : 1/0 ) with points ps 1 t "smooth", \
    filename using 5:( strcol(2) eq "smooth_inner" ? $4 : 1/0 ) with points ps 1 t "smooth inner", \
    filename using 5:( strcol(2) eq "smooth_wait" ? $4 : 1/0 ) with points ps 1 t "smooth wait", \
    filename using 5:( strcol(2) eq "smooth_outer" ? $4 : 1/0 ) with points ps 1 t "smooth outer", \
    filename using 5:( strcol(2) eq "smooth_res" ? $4 : 1/0 ) with points ps 1 t "smooth res", \
    filename using 5:( strcol(2) eq "smooth_res_inner" ? $4 : 1/0 ) with points ps 1 t "smooth res inner", \
    filename using 5:( strcol(2) eq "smooth_res_wait" ? $4 : 1/0 ) with points ps 1 t "smooth res wait", \
    filename using 5:( strcol(2) eq "smooth_res_col_bc" ? $4 : 1/0 ) with points ps 1 t "smooth res col bc", \
    filename using 5:( strcol(2) eq "smooth_res_outer" ? $4 : 1/0 ) with points ps 1 t "smooth res outer", \
    filename using 5:( strcol(2) eq "smooth_res_wait_set" ? $4 : 1/0 ) with points ps 1 t "smooth res wait set", \
    filename using 5:( strcol(2) eq "smooth_final" ? $4 : 1/0 ) with points ps 1 t "smooth final"

filename_out=filename."_seq.png"
set output filename_out
plot \
    filename using 5:( strcol(2) eq "smooth" ? $4 : 1/0 ) with points pt 1 ps 1 t "smooth", \
    filename using 5:( strcol(2) eq "smooth_inner" ? $4 : 1/0 ) with points pt 2 ps 1 t "smooth inner", \
    filename using 5:( strcol(2) eq "smooth_wait" ? $4 : 1/0 ) with points pt 3 ps 1 t "smooth wait", \
    filename using 5:( strcol(2) eq "smooth_outer" ? $4 : 1/0 ) with points pt 4 ps 1 t "smooth outer", \
    filename using 5:( strcol(2) eq "smooth_res" ? $4 : 1/0 ) with points pt 1 ps 1 t "smooth res" , \
    filename using 5:( strcol(2) eq "smooth_res_inner" ? $4 : 1/0 ) with points pt 2 ps 1 t "smooth res inner", \
    filename using 5:( strcol(2) eq "smooth_res_wait" ? $4 : 1/0 ) with points pt 3 ps 1 t "smooth res wait", \
    filename using 5:( strcol(2) eq "smooth_res_outer" ? $4 : 1/0 ) with points pt 4 ps 1 t "smooth res outer"

set ylabel "Flop/s"
filename_out=filename."_flops.png"
set output filename_out
plot \
    filename using 5:( strcol(2) eq "smooth" ? $6/$4 : 1/0 ) with points pt 1 ps 1 t "smooth", \
    filename using 5:( strcol(2) eq "smooth_inner" ? $6/$4 : 1/0 ) with points pt 2  ps 1 t "smooth inner", \
    filename using 5:( strcol(2) eq "smooth_outer" ? $6/$4 : 1/0 ) with points pt 3 ps 1 t "smooth outer", \
    filename using 5:( strcol(2) eq "smooth_res" ? $6/$4 : 1/0 ) with points pt 1 ps 1 t "smooth res", \
    filename using 5:( strcol(2) eq "smooth_res_inner" ? $6/$4 : 1/0 ) with points pt 2 ps 1 t "smooth res inner", \
    filename using 5:( strcol(2) eq "smooth_res_outer" ? $6/$4 : 1/0 ) with points pt 3 ps 1 t "smooth res outer"



