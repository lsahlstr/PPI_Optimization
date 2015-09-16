set xlabel "iRMSD (A)"
set ylabel "Interaction energy (kcal/mol)"

plot 'sys'.sys.'.dat' using 3:5 lt -1 lc rgb "black" linewidth 1.5 ps 1 pt 65

set lmargin 1
set size 1,1
set origin 0,0

#set yrange [-25:2]
#set xrange [1:50]

set key off
set size square
set origin 0,0

set terminal postscript eps color lw 2 "Helvetica" 30
set out 'sys'.sys.'_old.eps'
replot
set encoding iso_8859_1
set term pop
