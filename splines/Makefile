# Compiler and tools
CSC = mcs
GNUPLOT = gnuplot
SED = sed

all: spline_plot_A.png spline_plot_B.png

# Exercise A
ExerciseA.exe: ExerciseA.cs
	$(CSC) -out:$@ $<

spline_data_A.txt plot_spline_A.gp: ExerciseA.exe
	mono ExerciseA.exe
	mv spline_data.txt spline_data_A.txt
	mv plot_spline.gp plot_spline_A.gp
	$(SED) -i 's/spline_data.txt/spline_data_A.txt/g' plot_spline_A.gp

spline_plot_A.png: spline_data_A.txt plot_spline_A.gp
	$(GNUPLOT) plot_spline_A.gp
	mv spline_plot.png spline_plot_A.png

# Exercise B
ExerciseB.exe: ExerciseB.cs
	$(CSC) -out:$@ $<

spline_data_B.txt plot_spline_B.gp: ExerciseB.exe
	mono ExerciseB.exe
	mv spline_data.txt spline_data_B.txt
	mv plot_spline.gp plot_spline_B.gp
	$(SED) -i 's/spline_data.txt/spline_data_B.txt/g' plot_spline_B.gp

spline_plot_B.png: spline_data_B.txt plot_spline_B.gp
	$(GNUPLOT) plot_spline_B.gp
	mv spline_plot.png spline_plot_B.png

clean:
	rm -f *.exe spline_data_*.txt plot_spline_*.gp spline_plot_*.png

.PHONY: all clean
