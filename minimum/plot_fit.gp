set terminal svg size 800,600 enhanced
set output 'fit_plot.svg'
set title "Higgs data and fit"
set xlabel "Energy [GeV]"
set ylabel "Signal"
set grid

set style line 1 lc rgb "#8C3394" pt 7 ps 1.3 lw 1.5
set style line 2 lc rgb "#33A02C" lt 1 lw 2

plot \
  "exp_data.txt"  using 1:2:3 with yerrorbars ls 1 title 'Experimental Data', \
  "fit_curve.txt" using 1:2      with lines       ls 2 title 'Fitted Curve'
