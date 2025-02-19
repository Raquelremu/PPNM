CXXFLAGS = -Wall -std=c++23 -O3 -ffast-math
CODE = $(filter %.cs,$^)
LIBS = $(addprefix -reference:,$(filter %.dll,$^))
MKEXE = mcs -target:exe -out:$@ $(LIBS) $(CODE)
MKLIB = mcs -target:library -out:$@ $(LIBS) $(CODE)
TIME = time --output=$@ --append --format "$$nthreads %e %U"

times.svg : times.txt makefile
	echo '\
		set terminal svg background "white" ;\
		set output "$@" ;\
		set xlabel "number of threads" ;\
		set ylabel "running time" ;\
		set title "running time as function of number of threads" ;\
		plot "$<" with lines ;\
	' | tee log.gpi | gnuplot

times.txt: multiproc.exe #makefile
	>$@
	for nthreads in 1 2 3 4 5 6 7 8 9 10; \
		do $(TIME) mono multiproc.exe -threads $$nthreads -terms 2e9; \
		done

multiproc.exe: multiproc.cs ; $(MKEXE)
