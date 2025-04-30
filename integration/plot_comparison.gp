set terminal pngcairo enhanced font 'Arial,12'
set output 'comparison.png'
set style data histograms
set style fill solid
set boxwidth 0.8
set title 'Integration Results (with error estimates)'
set ylabel 'Value'
set xtics rotate by -45
set yrange [-5:3]
plot 'comparison_results.dat' using 2:xtic(1) title 'C#', 'comparison_results.dat' using 3 title 'Exact', '' using 0:($2>$3?$2:$3):4 with labels offset 0,1 font ',10' notitle