# Makefile for running ExerciseA.cs and plotting results

CS_COMPILER = mcs
EXECUTABLE = ExerciseA.exe
SOURCE = ExerciseA.cs
PLOT_SCRIPT = plot.gp
PLOTS = plot1.png plot2.png

all: run plot

build:
	$(CS_COMPILER) $(SOURCE)

run: build
	mono $(EXECUTABLE)

plot:
	gnuplot $(PLOT_SCRIPT)

clean:
	rm -f $(EXECUTABLE) results1.dat results2.dat $(PLOT_SCRIPT) $(PLOTS)

.PHONY: all build run plot clean

