
set terminal pngcairo enhanced font 'Arial,12'
set output 'circular_orbit.png'
set title 'Newtonian Circular Orbit'
set xlabel '{/Symbol f}'
set ylabel 'u({/Symbol f})'
set grid
plot 'circular_orbit.data' using 1:2 with lines title 'u({/Symbol f})'
