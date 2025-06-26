
set terminal pngcairo enhanced font 'Arial,12'
set output 'relativistic_orbit.png'
set title 'Relativistic Precession'
set xlabel '{/Symbol f}'
set ylabel 'u({/Symbol f})'
set grid
plot 'relativistic_orbit.data' using 1:2 with lines title 'u({/Symbol f})'

set output 'relativistic_orbit_polar.png'
set polar
set title 'Relativistic Precession (Polar Plot)'
unset xtics
unset ytics
set grid polar
set size square
plot 'relativistic_orbit.data' using 1:(1/$2) with lines title 'Orbit'
