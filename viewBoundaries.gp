set size square
set xrange [0:10]
set yrange [0:10]
set grid

set style line 1 \
    pointtype 7 pointsize 1.5 \
    linecolor rgb '#008000'
set style line 2 \
    linecolor rgb '#FF0000'

#nwall = 25 
#plot for[i = 1:nwall] 'wall'.i.'.txt' notitle with points linestyle 1, \
#     for[i = 1:nwall] 'wall'.i.'.txt' notitle with vectors linestyle 2

nbndr = 4
plot for[i = 1:nbndr] 'bndr'.i.'.txt' notitle with points linestyle 1, \
     for[i = 1:nbndr] 'bndr'.i.'.txt' notitle with vectors linestyle 2
