
set terminal pngcairo enhanced font 'Arial,12'
set output 'oscillator.png'
set title 'Simple Harmonic Oscillator'
set xlabel 'Time'
set ylabel 'Position and Velocity'
set grid
plot 'oscillator.data' using 1:2 with lines title 'Position (u)', \
     'oscillator.data' using 1:3 with lines title "Velocity (u')"
