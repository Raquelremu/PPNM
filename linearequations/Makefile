CFLAGS = -out:main.exe

all:
	mcs $(CFLAGS) QR.cs

run:
	mono main.exe

clean:
	rm -f main.exe out.times.data out err plot.png

out.times.data: main.exe
 >$@
 for N in $$(seq 100 20 200); do \
 time --format "$$N %e" --output $@ --append \
 mono $< -size:$$N 1>out 2>err ;\
 done 

plot: out.times.data
	echo '\
		set terminal png;\
		set output "Time.png";\
		set xlabel "N";\
		set ylabel "Time";\
		a=1; b=3 ;\
		f(x)=a*x**b ;\
		fit f(x) "$<" via a,b ;\
		plot "$<" with points \
		,f(x) with lines \
		' | gnuplot

#			,f(x) \
