
set terminal pngcairo enhanced font 'Arial,12'
set output 'elliptical_orbit.png'
set title 'Newtonian Elliptical Orbit'
set xlabel '{/Symbol f}'
set ylabel 'u({/Symbol f})'
set grid
plot 'elliptical_orbit.data' using 1:2 with lines title 'u({/Symbol f})'
