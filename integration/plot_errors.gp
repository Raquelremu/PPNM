set terminal pngcairo enhanced font 'Arial,12'
set output 'errors.png'
set logscale xy
set format x '10^{%L}'
set format y '10^{%L}'
set title 'Error Analysis (log-log scale)'
set xlabel 'Requested Accuracy'
set ylabel 'Error Magnitude'
set key top left
plot 'error_data.dat' using 1:2 with linespoints pt 7 title 'Estimated', 'error_data.dat' using 1:3 with linespoints pt 9 title 'Actual'