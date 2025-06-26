
set terminal pngcairo enhanced font 'arial,10' size 800,600
set output 'spline_plot.png'

set multiplot layout 2,1 title 'Quadratic Spline Interpolation'

set title 'Interpolation'
plot 'spline_data_B.txt' using 1:2 with lines title 'Qspline', \
     '' using 1:4 with lines title 'Exact sin(x)'

set title 'Integral'
plot 'spline_data_B.txt' using 1:3 with lines title 'Qspline Integral', \
     '' using 1:5 with lines title 'Exact integral (1-cos(x))'

unset multiplot
