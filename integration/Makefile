# Makefile for Numerical Integration Project

# Compiler and tools
CSC = mcs
GNUPLOT = gnuplot
RM = rm -f

# Targets
TARGET = integration.exe
SOURCE = integration.cs
PLOTS = errors.png evaluations.png comparison.png
DATA_FILES = error_data.dat comparison_results.dat
GP_SCRIPTS = plot_errors.gp plot_evaluations.gp plot_comparison.gp

# Default target
all: compile run plots

# Compilation
compile: $(TARGET)

$(TARGET): $(SOURCE)
	$(CSC) $(SOURCE) -out:$(TARGET)

# Execution
run: $(TARGET)
	mono $(TARGET)

# Plot generation
plots: $(PLOTS)

errors.png: plot_errors.gp error_data.dat
	$(GNUPLOT) plot_errors.gp

evaluations.png: plot_evaluations.gp error_data.dat
	$(GNUPLOT) plot_evaluations.gp

comparison.png: plot_comparison.gp comparison_results.dat
	$(GNUPLOT) plot_comparison.gp

# Cleanup
clean:
	$(RM) $(TARGET) $(PLOTS) $(DATA_FILES) $(GP_SCRIPTS)

.PHONY: all compile run plots clean
