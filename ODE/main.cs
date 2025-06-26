using System;
using static System.Math;
using static System.Console;
using System.IO;

public class tests {
    public static void Main(string[] args) {
   
        TestSimpleOscillator();
        TestExactSolution();
        
        TestNewtonianCircularOrbit();
        TestNewtonianEllipticalOrbit();
        TestRelativisticPrecession();
        
        GenerateGnuplotScripts();
        WriteLine("All tests completed. Run 'make plots' to generate the plots.");
    }


    public static void TestSimpleOscillator() {
        WriteLine("Testing simple harmonic oscillator u'' = -u...");

        Func<double, ode.vector, ode.vector> f = (x, y) => {
            return new ode.vector(y[1], -y[0]);
        };
        
        ode.vector y0 = new ode.vector(0.0, 1.0);
        var (xs, ys) = ode.driver(f, (0, 10), y0, h:0.1, acc:1e-6, eps:1e-6);
        
        using (StreamWriter writer = new StreamWriter("oscillator.data")) {
            for (int i = 0; i < xs.size; i++) {
                writer.WriteLine($"{xs[i]} {ys[i][0]} {ys[i][1]}");
            }
        }
        WriteLine("Saved oscillator solution to oscillator.data");
    }

    public static void TestExactSolution() {
        WriteLine("Testing exact solution for y'' = 2x...");
        
        Func<double, ode.vector, ode.vector> f = (x, y) => {
            return new ode.vector(y[1], 2*x);
        };
        
        ode.vector y0 = new ode.vector(0.0, 0.0);
        var (xs, ys) = ode.driver(f, (0, 10), y0, h:0.1, acc:1e-6, eps:1e-6, method:"rk23");
        
        using (StreamWriter writer = new StreamWriter("exact_solution.data")) {
            for (int i = 0; i < xs.size; i++) {
                double exact = xs[i]*xs[i]*xs[i]/3.0;
                writer.WriteLine($"{xs[i]} {ys[i][0]} {exact}");
            }
        }
        WriteLine("Saved exact solution test to exact_solution.data");
    }

    public static void TestNewtonianCircularOrbit() {
        WriteLine("Testing Newtonian circular orbit...");
        
        double ε = 0.0;
        Func<double, ode.vector, ode.vector> f = (φ, y) => {
            return new ode.vector(y[1], 1 - y[0] + ε*y[0]*y[0]);
        };
        
        ode.vector y0 = new ode.vector(1.0, 0.0); // u(0)=1, u'(0)=0
        var (φs, ys) = ode.driver(f, (0, 10*2*PI), y0, h:0.1, acc:1e-6, eps:1e-6);

        using (StreamWriter writer = new StreamWriter("circular_orbit.data")) {
            for (int i = 0; i < φs.size; i++) {
                writer.WriteLine($"{φs[i]} {ys[i][0]}");
            }
        }
        WriteLine("Saved circular orbit to circular_orbit.data");
    }


 public static void TestNewtonianEllipticalOrbit() {
    WriteLine("Testing Newtonian elliptical orbit...");
    
    const double ε = 0.0;
    const double initial_u_prime = -0.5;
    
    Func<double, ode.vector, ode.vector> f = (φ, y) => {
        return new ode.vector(y[1], 1 - y[0] + ε*y[0]*y[0]);
    };

    ode.vector y0 = new ode.vector(1.0, initial_u_prime);
    
    var (φs, ys) = ode.driver(
        f, 
        interval: (0, 6*2*PI),
        yinit: y0,
        h: 0.05,              
        acc: 1e-4,            
        eps: 1e-4,
        method: "rk12"        
    );
    
    using (StreamWriter writer = new StreamWriter("elliptical_orbit.data")) {
        for (int i = 0; i < φs.size; i++) {
            writer.WriteLine($"{φs[i]} {ys[i][0]}");
        }
    }
    WriteLine("Saved elliptical orbit to elliptical_orbit.data");
}

    public static void TestRelativisticPrecession() {
        WriteLine("Testing relativistic precession...");
        
        double ε = 0.01;
        Func<double, ode.vector, ode.vector> f = (φ, y) => {
            return new ode.vector(y[1], 1 - y[0] + ε*y[0]*y[0]);
        };
        
        ode.vector y0 = new ode.vector(1.0, -0.5); 
        var (φs, ys) = ode.driver(f, (0, 20*2*PI), y0, h:0.1, acc:1e-6, eps:0.01);

        using (StreamWriter writer = new StreamWriter("relativistic_orbit.data")) {
            for (int i = 0; i < φs.size; i++) {
                writer.WriteLine($"{φs[i]} {ys[i][0]}");
            }
        }
        WriteLine("Saved relativistic orbit to relativistic_orbit.data");
    }
 

    public static void GenerateGnuplotScripts() {
        File.WriteAllText("plot_oscillator.gp", @"
set terminal pngcairo enhanced font 'Arial,12'
set output 'oscillator.png'
set title 'Simple Harmonic Oscillator'
set xlabel 'Time'
set ylabel 'Position and Velocity'
set grid
plot 'oscillator.data' using 1:2 with lines title 'Position (u)', \
     'oscillator.data' using 1:3 with lines title ""Velocity (u')""
");

        File.WriteAllText("plot_exact_solution.gp", @"
set terminal pngcairo enhanced font 'Arial,12'
set output 'exact_solution.png'
set title ""Exact Solution Test (y'' = 2x)""
set xlabel 'x'
set ylabel 'y(x)'
set grid
plot 'exact_solution.data' using 1:2 with lines title 'Numerical', \
     'exact_solution.data' using 1:3 with lines title 'Exact (x^3/3)'
");

        File.WriteAllText("plot_circular_orbit.gp", @"
set terminal pngcairo enhanced font 'Arial,12'
set output 'circular_orbit.png'
set title 'Newtonian Circular Orbit'
set xlabel '{/Symbol f}'
set ylabel 'u({/Symbol f})'
set grid
plot 'circular_orbit.data' using 1:2 with lines title 'u({/Symbol f})'
");

        File.WriteAllText("plot_elliptical_orbit.gp", @"
set terminal pngcairo enhanced font 'Arial,12'
set output 'elliptical_orbit.png'
set title 'Newtonian Elliptical Orbit'
set xlabel '{/Symbol f}'
set ylabel 'u({/Symbol f})'
set grid
plot 'elliptical_orbit.data' using 1:2 with lines title 'u({/Symbol f})'
");

        File.WriteAllText("plot_relativistic_orbit.gp", @"
set terminal pngcairo enhanced font 'Arial,12'
set output 'relativistic_orbit.png'
set title 'Relativistic Precession'
set xlabel '{/Symbol f}'
set ylabel 'u({/Symbol f})'
set grid
plot 'relativistic_orbit.data' using 1:2 with lines title 'u({/Symbol f})'

set output 'relativistic_orbit_polar.png'
set polar
set title 'Relativistic Precession (Polar Plot)'
unset xtics
unset ytics
set grid polar
set size square
plot 'relativistic_orbit.data' using 1:(1/$2) with lines title 'Orbit'
");
    }
}
