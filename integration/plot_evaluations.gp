set terminal pngcairo enhanced font 'Arial,12'
set output 'evaluations.png'
set logscale xy
set format x '10^{%L}'
set format y '%g'
set title 'Function Evaluations vs Accuracy'
set xlabel 'Requested Accuracy'
set ylabel 'Evaluations'
plot 'error_data.dat' using 1:4 with linespoints pt 5 lw 2 title ''