# Gnuplot script file for plotting data 
# This file is called   plot.p
set terminal postscript eps enhanced color
set output "res2.eps"
set key Left reverse
set style line 1 lc rgb '#0000ff' pt 6 ps 1.5 lt 1 lw 5
set style line 2 lc rgb '#8b1a0e' pt 1 ps 1.5 lt 1 lw 5
set style line 3 lc rgb '#000000' pt 10 ps 1.5 lt 1 lw 5
set xrange [0:300]
set yrange [60:150]
#set ytics nomirror
#set y2tics
set xlabel "CPU time"
set ylabel "Solution Cost"
   plot './outfile/data-abeast2' u 1:2 w lp ls 1 t 'abeast' axes x1y1, \
   './outfile/data-abeastcost2' u 1:2 w lp ls 2 t 'abeast-onlyPruning' axes x1y1, './outfile/data-sst2' u 1:2 w lp ls 3 t 'sststar' axes x1y1
set output


