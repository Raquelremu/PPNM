
set terminal pngcairo enhanced font 'arial,10' size 800,600
set output 'spline_plot.png'

set multiplot layout 2,1 title 'Linear Spline Interpolation'

set title 'Interpolation'
plot 'spline_data_A.txt' using 1:2 with lines title 'Linear spline', \
     '' using 1:4 with lines title 'Exact cos(x)'

set title 'Integral'
plot 'spline_data_A.txt' using 1:3 with lines title 'Spline integral', \
     '' using 1:5 with lines title 'Exact integral (sin)'

unset multiplot
