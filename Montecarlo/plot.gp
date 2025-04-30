
set terminal pngcairo enhanced font 'arial,10' size 1200,800

# First plot: Unit circle results
set output 'plot1.png'
set multiplot layout 2,1

set logscale xy
set xlabel 'Number of points (N)'
set ylabel 'Error'
set title 'Monte Carlo Integration - Unit Circle Area'
plot 'results1.dat' using 1:3 with linespoints title 'Estimated Error', \
     'results1.dat' using 1:4 with linespoints title 'Actual Error', \
     1/sqrt(x) with lines title '1/sqrt(N) scaling'

unset logscale
set ylabel 'Estimated π value'
set title 'Convergence to π'
plot 'results1.dat' using 1:2 with linespoints title 'Monte Carlo', \
     pi with lines title 'Exact π value'

unset multiplot

# Second plot: Singular integral results
set output 'plot2.png'
set multiplot layout 2,1

set logscale xy
set xlabel 'Number of points (N)'
set ylabel 'Error'
set title 'Monte Carlo Integration - Singular Integral'
plot 'results2.dat' using 1:3 with linespoints title 'Estimated Error', \
     'results2.dat' using 1:4 with linespoints title 'Actual Error'

unset logscale
set ylabel 'Estimated value'
set title 'Convergence to exact value'
exact_value = 1.3932039296856768591842462603255
plot 'results2.dat' using 1:2 with linespoints title 'Monte Carlo', \
     exact_value with lines title 'Exact value'

unset multiplot
