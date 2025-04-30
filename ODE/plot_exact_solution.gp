
set terminal pngcairo enhanced font 'Arial,12'
set output 'exact_solution.png'
set title "Exact Solution Test (y'' = 2x)"
set xlabel 'x'
set ylabel 'y(x)'
set grid
plot 'exact_solution.data' using 1:2 with lines title 'Numerical', \
     'exact_solution.data' using 1:3 with lines title 'Exact (x^3/3)'
