CSC       := mcs
CSFILES   := $(wildcard *.cs)
TARGET    := assignment.exe
INPUT     := input.csv
PLOT      := evaluations.png

.PHONY: all run plot clean

all: $(TARGET)

$(TARGET): $(CSFILES)
	$(CSC) -debug -out:$(TARGET) $(CSFILES)

# generate sample input if missing
$(INPUT):
	@echo "x,y" >  $(INPUT)
	@echo "0,0" >> $(INPUT)
	@echo "1,1" >> $(INPUT)
	@echo "2,0.5" >> $(INPUT)
	@echo "3,1.5" >> $(INPUT)
	@echo "4,1" >> $(INPUT)
	@echo "5,2" >> $(INPUT)
	@echo "Generated sample $(INPUT)"

# run the program
run: all $(INPUT)
	@mono $(TARGET)
	@$(MAKE) plot

# plot target
plot: evaluations.csv
	@echo "Plotting evaluations.csv → $(PLOT)"
	@gnuplot -e "\
	set datafile separator ','; \
	set terminal png size 800,600; \
	set output '$(PLOT)'; \
	set title 'Akima Spline Interpolation'; \
	set xlabel 'x'; set ylabel 'Interpolated y'; \
	plot 'evaluations.csv' using 1:2 with linespoints lw 2 pt 7 title 'y(z)'\
	"
	@echo "Saved plot as $(PLOT)"

clean:
	@rm -f $(TARGET) input_data.csv evaluations.csv plot.gp $(PLOT)

