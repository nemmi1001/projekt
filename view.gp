set xrange [0:10]
set yrange [0:10]
set grid

set style line 1 \
    pointtype 7 pointsize 1 \
    linecolor rgb '#008000'
set style line 2 \
    linecolor rgb '#FF0000'

#plot "wall.txt" notitle with linespoints linestyle 1, \
#     "wall.txt" notitle with vectors linestyle 2

plot "wall1.txt" notitle with linespoints linestyle 1, \
     "wall1.txt" notitle with vectors linestyle 2, \
     "wall2.txt" notitle with linespoints linestyle 1, \
     "wall2.txt" notitle with vectors linestyle 2, \
     "wall3.txt" notitle with linespoints linestyle 1, \
     "wall3.txt" notitle with vectors linestyle 2, \
     "wall4.txt" notitle with linespoints linestyle 1, \
     "wall4.txt" notitle with vectors linestyle 2