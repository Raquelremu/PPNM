using System;
using System.IO;

public class LinearSpline
{
    public static int BinSearch(double[] x, double z)
    {
        if (z < x[0] || z > x[x.Length - 1])
            throw new ArgumentException("binsearch: z is outside the x range");
        
        int i = 0, j = x.Length - 1;
        while (j - i > 1)
        {
            int mid = (i + j) / 2;
            if (z > x[mid]) i = mid;
            else j = mid;
        }
        return i;
    }

    public static double Linterp(double[] x, double[] y, double z)
    {
        if (x.Length != y.Length)
            throw new ArgumentException("x and y arrays must have the same length");
        
        int i = BinSearch(x, z);
        double dx = x[i + 1] - x[i];
        if (!(dx > 0)) throw new ArgumentException("x values must be strictly increasing");
        
        double dy = y[i + 1] - y[i];
        return y[i] + dy / dx * (z - x[i]);
    }

    public static double LinterpIntegral(double[] x, double[] y, double z)
    {
        if (x.Length != y.Length)
            throw new ArgumentException("x and y arrays must have the same length");
        
        if (z < x[0] || z > x[x.Length - 1])
            throw new ArgumentException("z is outside the x range");

        int i = BinSearch(x, z);
        double integral = 0.0;

        for (int j = 0; j < i; j++)
        {
            double dx = x[j + 1] - x[j];
            double dy = y[j + 1] - y[j];
            double p = dy / dx;
            integral += y[j] * dx + p * dx * dx / 2;
        }

        integral += y[i] * (z - x[i]) + ((Linterp(x, y, z) - y[i])/(z - x[i])) * (z - x[i]) * (z - x[i]) / 2;

        return integral;
    }
}

class Program
{
    static void Main()
    {
        int n = 10;
        double[] x = new double[n];
        double[] y = new double*[n];
        
        for (int i = 0; i < n; i++)
        {
            x[i] = i;
            y[i] = Math.Cos(i);
        }

        using (StreamWriter sw = new StreamWriter("spline_data.txt"))
        {
            sw.WriteLine("# z\tInterpolated\tIntegral\tExactCos\tExactIntegral");
            for (double z = x[0]; z <= x[n-1]; z += 0.1)
            {
                double interpolated = LinearSpline.Linterp(x, y, z);
                double integral = LinearSpline.LinterpIntegral(x, y, z);
                sw.WriteLine($"{z}\t{interpolated}\t{integral}\t{Math.Cos(z)}\t{Math.Sin(z) - Math.Sin(x[0])}");
            }
        }

        File.WriteAllText("plot_spline.gp", @"
set terminal pngcairo enhanced font 'arial,10' size 800,600
set output 'spline_plot.png'

set multiplot layout 2,1 title 'Linear Spline Interpolation'

set title 'Interpolation'
plot 'spline_data.txt' using 1:2 with lines title 'Linear spline', \
     '' using 1:4 with lines title 'Exact cos(x)'

set title 'Integral'
plot 'spline_data.txt' using 1:3 with lines title 'Spline integral', \
     '' using 1:5 with lines title 'Exact integral (sin)'

unset multiplot
");
    }
}
