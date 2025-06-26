# Makefile that compiles and runs ExerciseA.cs,
# then splits its output into two files (exp_data.txt and fit_data.txt) for plotting

# Variables
EXE = ExerciseA.exe
SRC = ExerciseA.cs
DATA = plot_data.txt
PNG = plot.png

.PHONY: all clean

all: $(PNG)

# Compile the C# source file using mcs (Mono C# compiler)
$(EXE): $(SRC)
	mcs $(SRC) -out:$(EXE)

# Run the executable to produce plot_data.txt
$(DATA): $(EXE)
	mono $(EXE)

# Create the plot by splitting plot_data.txt into two separate files for the two blocks,
# and then using a temporary gnuplot script.
$(PNG): $(DATA)
	@echo "Extracting experimental data..."
	@awk 'BEGIN {p=1} /^$$/ {p=0} {if(p && $$0 !~ /^#/ && NF > 0) print}' $(DATA) > exp_data.txt
	@echo "Extracting fitted data..."
	@awk 'BEGIN {p=0} /^$$/ {p=1; next} {if(p && $$0 !~ /^#/ && NF > 0) print}' $(DATA) > fit_data.txt
	@echo "Creating gnuplot script..."
	@echo "set terminal pngcairo size 800,600" > plot.gp
	@echo "set output '$(PNG)'" >> plot.gp
	@echo "set title 'Experimental Data and Fitted Curve'" >> plot.gp
	@echo "set xlabel 'Time (t)'" >> plot.gp
	@echo "set ylabel 'Activity (y)'" >> plot.gp
	@echo "set key left top" >> plot.gp
	@echo "plot 'exp_data.txt' using 1:2:3 with yerrorbars title 'Experimental Data', \\" >> plot.gp
	@echo "     'fit_data.txt' using 1:2 with lines title 'Fitted Curve'" >> plot.gp
	@echo "Running gnuplot..."
	gnuplot plot.gp
	@rm plot.gp exp_data.txt fit_data.txt
	@echo "Plot generated: $(PNG)"

# Clean up generated files
clean:
	rm -f $(EXE) $(DATA) $(PNG)

