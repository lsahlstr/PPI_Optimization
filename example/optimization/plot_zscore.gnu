set xlabel "cycle"
set ylabel "Z-score"

plot "Zscore.txt" using 1:2 linetype 1 lw 4 ps 1 pt 6 lc rgb "black" with lines

set lmargin 1
set size 1,1
set origin 0,0

#set yrange [-25:2]
#set xrange [1:50]
set xtics 0,30000,100000

set key off
set size square
set origin 0,0

set terminal postscript eps color lw 2 "Helvetica" 24
set out 'zscore.eps'
replot
set encoding iso_8859_1
set term pop
