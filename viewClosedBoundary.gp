set size square
set xrange [0:10]
set yrange [0:10]
set grid

set style line 1 \
    pointtype 7 pointsize 1.5 \
    linecolor rgb '#008000'
set style line 2 \
    linecolor rgb '#FF0000'

plot "closed_boundary.txt" notitle with points linestyle 1, \
     "closed_boundary.txt" notitle with vectors linestyle 2 