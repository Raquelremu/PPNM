set terminal pngcairo enhanced font 'Arial,14'
set output 'plot.png'
set xlabel 'Matrix Size (N)'
set ylabel 'Execution Time (seconds)'
set title 'QR Decomposition Time Complexity'
set grid
f(x) = a*x**3 + b
fit f(x) 'out.times.data' using 1:2 via a, b
plot 'out.times.data' using 1:2 with points title 'Measured Time', \
     f(x) with lines title 'Fitted O(N^3)'
