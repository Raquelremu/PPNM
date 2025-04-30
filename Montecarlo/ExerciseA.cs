using System;
using System.Linq;
using System.IO;
using System.Diagnostics;

public class vector
{
    private double[] data;

    public vector(int size) => data = new double[size];
    public vector(double[] data) => this.data = data;

    public int Size => data.Length;
    public double this[int i]
    {
        get => data[i];
        set => data[i] = value;
    }

    public static implicit operator vector(double[] array) => new vector(array);
}

public class MonteCarlo
{
    public static Tuple<double, double> PlainMC(Func<vector, double> f, vector a, vector b, int N) {
        int dim = a.Size;
        double V = 1;
        for (int i = 0; i < dim; i++) V *= b[i] - a[i];

        double sum = 0, sum2 = 0;
        vector x = new vector(dim);
        Random rnd = new Random();

        for (int i = 0; i < N; i++) {
            for (int k = 0; k < dim; k++)
                x[k] = a[k] + rnd.NextDouble() * (b[k] - a[k]);
            double fx = f(x);
            sum += fx;
            sum2 += fx * fx;
        }

        double mean = sum / N;
        double sigma = Math.Sqrt(sum2 / N - mean * mean);
        return Tuple.Create(mean * V, sigma * V / Math.Sqrt(N));
    }
}

public class Program
{
    public static void Main()
    {
        Console.WriteLine("Running unit circle area test...");
        RunUnitCircleTest();

        Console.WriteLine("\nRunning singular integral test...");
        RunSingularIntegralTest();

        GeneratePlots();
    }

    private static void RunUnitCircleTest()
    {
        var a = new vector(new[] { -1.0, -1.0 });
        var b = new vector(new[] { 1.0, 1.0 });

        // Function for computing the unit circle area: points inside the circle contribute 1; outside contribute 0.
        Func<vector, double> circleFunc = v =>
            (v[0] * v[0] + v[1] * v[1] <= 1) ? 1.0 : 0.0;

        int[] Ns = { 100, 1000, 10000, 100000, 1000000 };
        var results = Ns.Select(N =>
        {
            var tuple = MonteCarlo.PlainMC(circleFunc, a, b, N);
            double result = tuple.Item1;
            double error = tuple.Item2;
            double actualError = Math.Abs(result - Math.PI);
            return new Tuple<int, double, double, double>(N, result, error, actualError);
        }).ToList();

        using (var writer = new StreamWriter("results1.dat"))
        {
            writer.WriteLine("N\tResult\tEstimatedError\tActualError");
            foreach (var r in results)
            {
                writer.WriteLine($"{r.Item1}\t{r.Item2}\t{r.Item3}\t{r.Item4}");
            }
        }

        Console.WriteLine("Unit circle test results saved to results1.dat");
    }

    private static void RunSingularIntegralTest()
    {
        var a = new vector(new[] { 0.0, 0.0, 0.0 });
        var b = new vector(new[] { Math.PI, Math.PI, Math.PI });

        // Improved integrand for singular integral with careful handling near singularities.
        Func<vector, double> singularFunc = v =>
        {
            double cosx = Math.Cos(v[0]);
            double cosy = Math.Cos(v[1]);
            double cosz = Math.Cos(v[2]);
            double product = cosx * cosy * cosz;

            // Skip points too close to the singular region.
            if (1 - product < 1e-12)
                return 0;

            return 1.0 / (1.0 - product);
        };

        double exactValue = 1.3932039296856768591842462603255;
        int[] Ns = { 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000 };

        var results = Ns.Select(N =>
        {
            var tuple = MonteCarlo.PlainMC(singularFunc, a, b, N);
            double result = tuple.Item1;
            double error = tuple.Item2;
            double actualError = Math.Abs(result - exactValue);
            return new Tuple<int, double, double, double>(N, result, error, actualError);
        }).ToList();

        using (var writer = new StreamWriter("results2.dat"))
        {
            writer.WriteLine("N\tResult\tEstimatedError\tActualError");
            foreach (var r in results)
            {
                writer.WriteLine($"{r.Item1}\t{r.Item2}\t{r.Item3}\t{r.Item4}");
            }
        }

        Console.WriteLine("Singular integral test results saved to results2.dat");
    }

    private static void GeneratePlots()
    {
        string gnuplotScript = @"
set terminal pngcairo enhanced font 'arial,10' size 1200,800

# First plot: Unit circle results
set output 'plot1.png'
set multiplot layout 2,1

set logscale xy
set xlabel 'Number of points (N)'
set ylabel 'Error'
set title 'Monte Carlo Integration - Unit Circle Area'
plot 'results1.dat' using 1:3 with linespoints title 'Estimated Error', \
     'results1.dat' using 1:4 with linespoints title 'Actual Error', \
     1/sqrt(x) with lines title '1/sqrt(N) scaling'

unset logscale
set ylabel 'Estimated π value'
set title 'Convergence to π'
plot 'results1.dat' using 1:2 with linespoints title 'Monte Carlo', \
     pi with lines title 'Exact π value'

unset multiplot

# Second plot: Singular integral results
set output 'plot2.png'
set multiplot layout 2,1

set logscale xy
set xlabel 'Number of points (N)'
set ylabel 'Error'
set title 'Monte Carlo Integration - Singular Integral'
plot 'results2.dat' using 1:3 with linespoints title 'Estimated Error', \
     'results2.dat' using 1:4 with linespoints title 'Actual Error'

unset logscale
set ylabel 'Estimated value'
set title 'Convergence to exact value'
exact_value = 1.3932039296856768591842462603255
plot 'results2.dat' using 1:2 with linespoints title 'Monte Carlo', \
     exact_value with lines title 'Exact value'

unset multiplot
";

        File.WriteAllText("plot.gp", gnuplotScript);

        var gnuplot = new Process
        {
            StartInfo = new ProcessStartInfo
            {
                FileName = "gnuplot",
                Arguments = "plot.gp",
                RedirectStandardOutput = true,
                UseShellExecute = false,
                CreateNoWindow = true,
            }
        };

        gnuplot.Start();
        gnuplot.WaitForExit();

        Console.WriteLine("\nPlots generated:");
        Console.WriteLine("plot1.png - Unit circle area results");
        Console.WriteLine("plot2.png - Singular integral results");
    }
}
