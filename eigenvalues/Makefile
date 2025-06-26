# Files and compiler settings
CSFILES   = Program.cs
TARGET    = main.exe
MONO      = mcs
CFLAGS    = -debug -out:$(TARGET)

all: $(TARGET)

$(TARGET): $(CSFILES)
	$(MONO) $(CFLAGS) $(CSFILES)


# Experiment 1
RMAX_FIXED = 10
DR_VALUES  = 0.1 0.15 0.2 0.25 0.3 0.35 0.4

# Experiment 2
DR_FIXED   = 0.3
RMAX_VALUES = 8 10 12 14 16

conv_dr: $(TARGET)
	@echo "Running convergence analysis: fixed rmax=$(RMAX_FIXED), varying dr"
	@rm -f out_dr.txt
	@for dr in $(DR_VALUES); do \
	  echo "Running for dr=$$dr"; \
	  time -ao out_dr.txt -f "$$dr %e" mono $(TARGET) -rmax $(RMAX_FIXED) -dr $$dr > conv_dr_$$dr.txt & \
	done; \
	wait; \
	sort out_dr.txt -o out_dr.txt
	@echo "Convergence (dr): timing information is in out_dr.txt"

conv_rmax: $(TARGET)
	@echo "Running convergence analysis: fixed dr=$(DR_FIXED), varying rmax"
	@rm -f out_rmax.txt
	@for rmax in $(RMAX_VALUES); do \
	  echo "Running for rmax=$$rmax"; \
	  time -ao out_rmax.txt -f "$$rmax %e" mono $(TARGET) -rmax $$rmax -dr $(DR_FIXED) > conv_rmax_$$rmax.txt & \
	done; \
	wait; \
	sort out_rmax.txt -o out_rmax.txt
	@echo "Convergence (rmax): timing information is in out_rmax.txt"

conv: conv_dr conv_rmax

plot_dr:
	@gnuplot -e "set terminal png; set output 'convergence_dr.png'; set title 'Convergence vs. Δr (rmax=$(RMAX_FIXED))'; set xlabel 'Δr'; set ylabel 'Time (sec)'; plot 'out_dr.txt' using 1:2 with linespoints title 'Δr'" 

plot_rmax:
	@gnuplot -e "set terminal png; set output 'convergence_rmax.png'; set title 'Convergence vs. rmax (Δr=$(DR_FIXED))'; set xlabel 'rmax'; set ylabel 'Time (sec)'; plot 'out_rmax.txt' using 1:2 with linespoints title 'rmax'"

clean:
	rm -f $(TARGET) out_dr.txt out_rmax.txt conv_dr_*.txt conv_rmax_*.txt *.png

.PHONY: all conv_dr conv_rmax conv plot_dr plot_rmax clean

