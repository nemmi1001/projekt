#!/bin/bash
make clean
make
./output
#//gnuplot
gnuplot  "view.gp" -p
