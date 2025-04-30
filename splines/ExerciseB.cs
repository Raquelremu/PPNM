using System;
using System.IO;

public class Qspline
{
    private double[] x, y, b, c;

    public Qspline(double[] xs, double[] ys)
    {
        if (xs.Length != ys.Length)
            throw new ArgumentException("x and y are not the same length");
        if (xs.Length < 2)
            throw new ArgumentException("Need at least 2 points");

        int n = xs.Length;
        x = new double[n];
        y = new double[n];
        Array.Copy(xs, x, n);
        Array.Copy(ys, y, n);

        b = new double[n - 1];
        c = new double[n - 1];

        double[] h = new double[n - 1];
        double[] p = new double[n - 1];

        for (int i = 0; i < n - 1; i++)
        {
            h[i] = x[i + 1] - x[i];
            p[i] = (y[i + 1] - y[i]) / h[i];
        }

        c[0] = 0;
        for (int i = 0; i < n - 2; i++)
        {
            c[i + 1] = (p[i + 1] - p[i] - c[i] * h[i]) / h[i + 1];
        }

        c[n - 2] /= 2;
        for (int i = n - 3; i >= 0; i--)
        {
            c[i] = (p[i + 1] - p[i] - c[i + 1] * h[i + 1]) / h[i];
        }

        for (int i = 0; i < n - 1; i++)
        {
            b[i] = p[i] - c[i] * h[i];
        }
    }

    public double Evaluate(double z)
    {
        int i = BinSearch(x, z);
        double h = z - x[i];
        return y[i] + h * (b[i] + h * c[i]);
    }

    public double Derivative(double z)
    {
        int i = BinSearch(x, z);
        double h = z - x[i];
        return b[i] + 2 * c[i] * h;
    }

    public double Integral(double z)
    {
        int i = BinSearch(x, z);
        double integral = 0.0;

        for (int j = 0; j < i; j++)
        {
            double h = x[j + 1] - x[j];
            integral += y[j] * h + b[j] * h * h / 2 + c[j] * h * h * h / 3;
        }

        double partialH = z - x[i];
        integral += y[i] * partialH + b[i] * partialH * partialH / 2 + c[i] * partialH * partialH * partialH / 3;
        
        return integral; 
    }

    public static int BinSearch(double[] x, double z)
    {
        if (z < x[0] || z > x[x.Length - 1])
        {
            throw new Exception("z is out of range");
        }
        int i = 0, j = x.Length - 1;
        while (j - i > 1)
        {
            int mid = (i + j) / 2;
            if (z > x[mid])
            {
                i = mid;
            }
            else
            {
                j = mid;
            }
        }
        return i;
    }
}

class Program
{
    public static int Main(string[] args)
    {
        int n = 10;
        double[] x = new double[n];
        double[] y = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = i;
            y[i] = Math.Sin(i);
        }
        Qspline s = new Qspline(x, y);

        using (StreamWriter sw = new StreamWriter("spline_data.txt"))
        {
            sw.WriteLine("# z\tSpline\tIntegral\tExactSin\tExactIntegral");
            for (double z = x[0]; z <= x[n - 1]; z += 0.1)
            {
                double splineVal = s.Evaluate(z);
                double integral = s.Integral(z);
                double exactSin = Math.Sin(z);
                double exactIntegral = 1 - Math.Cos(z);
                sw.WriteLine($"{z}\t{splineVal}\t{integral}\t{exactSin}\t{exactIntegral}");
            }
        }


        File.WriteAllText("plot_spline.gp", @"
set terminal pngcairo enhanced font 'arial,10' size 800,600
set output 'spline_plot.png'

set multiplot layout 2,1 title 'Quadratic Spline Interpolation'

set title 'Interpolation'
plot 'spline_data.txt' using 1:2 with lines title 'Qspline', \
     '' using 1:4 with lines title 'Exact sin(x)'

set title 'Integral'
plot 'spline_data.txt' using 1:3 with lines title 'Qspline Integral', \
     '' using 1:5 with lines title 'Exact integral (1-cos(x))'

unset multiplot
");

        return 0;
    }
}

