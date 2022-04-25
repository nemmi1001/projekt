set size square
set xrange [-5:10]
set yrange [-5:10]
set grid

set style line 1 \
    pointtype 7 pointsize 1.5 \
    linecolor rgb '#008000'
set style line 2 \
    linecolor rgb '#FF0000'

# points / lines / linespoints

nwall = 12
plot for[i = 1:nwall] 'wall'.i.'.txt' notitle with points linestyle 1, \
     for[i = 1:nwall] 'wall'.i.'.txt' notitle with vectors linestyle 2
