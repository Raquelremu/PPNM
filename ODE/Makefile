# Makefile for ODE homework

# Compiler settings
CSC = mcs
CSFLAGS = -target:library -out:ode.dll
TESTFLAGS = -target:exe -out:main.exe -reference:ode.dll

# Plotting command
PLOT = gnuplot

all: main.exe plots

ode.dll: ode.cs
	$(CSC) $(CSFLAGS) ode.cs

main.exe: main.cs ode.dll
	$(CSC) $(TESTFLAGS) main.cs

run: main.exe
	mono main.exe

plots: run
	@echo "Generating plots..."
	@$(PLOT) plot_oscillator.gp
	@$(PLOT) plot_exact_solution.gp
	@$(PLOT) plot_circular_orbit.gp
	@$(PLOT) plot_elliptical_orbit.gp
	@$(PLOT) plot_relativistic_orbit.gp
	@echo "All plots generated!"

clean:
	rm -f *.exe *.dll *.data *.png *.gp

.PHONY: all run plots clean
