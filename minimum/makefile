# Makefile

all: ExerciseAB.exe

ExerciseAB.exe: ExerciseAB.cs
	mcs ExerciseAB.cs -out:ExerciseAB.exe


plot_fit.gp:
	@echo "set terminal svg size 800,600 enhanced"                >  $@
	@echo "set output 'fit_plot.svg'"                             >> $@
	@echo "set title \"Higgs data and fit\""                      >> $@
	@echo "set xlabel \"Energy [GeV]\""                           >> $@
	@echo "set ylabel \"Signal\""                                 >> $@
	@echo "set grid"                                              >> $@
	@echo ""                                                     >> $@
	@echo "set style line 1 lc rgb \"#8C3394\" pt 7 ps 1.3 lw 1.5" >> $@
	@echo "set style line 2 lc rgb \"#33A02C\" lt 1 lw 2"         >> $@
	@echo ""                                                     >> $@
	@echo "plot \\"                                               >> $@
	@echo "  \"exp_data.txt\"  using 1:2:3 with yerrorbars ls 1 title 'Experimental Data', \\" >> $@
	@echo "  \"fit_curve.txt\" using 1:2      with lines       ls 2 title 'Fitted Curve'"       >> $@


run: ExerciseAB.exe

	mono ExerciseAB.exe < higgs.data.txt 1> out.txt 2> log
	

	mono ExerciseAB.exe < higgs.data.txt 2> log | tee out.txt
	
	$(MAKE) plot_fit.gp

	gnuplot plot_fit.gp
	display fit_plot.svg &


clean:
	rm -f ExerciseAB.exe out.txt log \
	      exp_data.txt fit_curve.txt \
	      plot_fit.gp fit_plot.svg

.PHONY: all run clean plot_fit.gp

