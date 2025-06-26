using System;
using System.Runtime.InteropServices;
using System.IO;
using System.Globalization;
using System.Collections.Generic;

public class AkimaSpline
{
    public int n;
    public double[] x;
    public double[] y;
    public double[] b;
    public double[] c;
    public double[] d;

    public AkimaSpline(int n, double[] x, double[] y)
    {
        System.Diagnostics.Debug.Assert(n > 2);
        this.n = n;
        this.x = new double[n];
        this.y = new double[n];
        this.b = new double[n];
        this.c = new double[n - 1];
        this.d = new double[n - 1];

        double[] h = new double[n - 1];
        double[] p = new double[n - 1];

        for (int i = 0; i < n - 1; i++)
        {
            h[i] = x[i + 1] - x[i];
            System.Diagnostics.Debug.Assert(h[i] > 0);
            p[i] = (y[i + 1] - y[i]) / h[i];
        }

        for (int i = 0; i < n; i++)
        {
            this.x[i] = x[i];
            this.y[i] = y[i];
        }

        this.b[0] = p[0];
        this.b[1] = (p[0] + p[1]) / 2;
        this.b[n - 1] = p[n - 2];
        this.b[n - 2] = (p[n - 2] + p[n - 3]) / 2;

        for (int i = 2; i < n - 2; i++)
        {
            double w1 = Math.Abs(p[i + 1] - p[i]);
            double w2 = Math.Abs(p[i - 1] - p[i - 2]);

            if (w1 + w2 == 0)
            {
                this.b[i] = (p[i - 1] + p[i]) / 2;
            }
            else
            {
                this.b[i] = (w1 * p[i - 1] + w2 * p[i]) / (w1 + w2);
            }
        }

        for (int i = 0; i < n - 1; i++)
        {
            this.c[i] = (3 * p[i] - 2 * this.b[i] - this.b[i + 1]) / h[i];
            this.d[i] = (this.b[i + 1] + this.b[i] - 2 * p[i]) / (h[i] * h[i]);
        }
    }

    public double Evaluate(double z)
    {
        System.Diagnostics.Debug.Assert(z >= this.x[0] && z <= this.x[this.n - 1]);

        int i = 0;
        int j = this.n - 1;
        
        while (j - i > 1)
        {
            int m = (i + j) / 2;
            if (z >= this.x[m])
                i = m;
            else
                j = m;
        }

        double h = z - this.x[i];
        return this.y[i] + h * (this.b[i] + h * (this.c[i] + h * this.d[i]));
    }
}

public class Program
{
    public static void Main()
    {
        var all = File.ReadAllLines("input.csv");
        File.WriteAllLines("input_data.csv", all);              
        int n = all.Length - 1;                                 
        double[] x = new double[n], y = new double[n];
        for (int i = 1; i < all.Length; i++)
        {
            var p = all[i].Split(',');
            x[i - 1] = double.Parse(p[0], CultureInfo.InvariantCulture);
            y[i - 1] = double.Parse(p[1], CultureInfo.InvariantCulture);
        }

        if (x.Length != y.Length)
        {
            Console.WriteLine("Error");
            Environment.Exit(1);
        }

        var spline = new AkimaSpline(n, x, y);

        Console.WriteLine("Akima Sub-spline Interpolation");
        Console.WriteLine("==============================");
        Console.WriteLine("x\tInterpolated y");
        Console.WriteLine("-----------------------");

        int m = 2 * n;
        var pts = new double[m + 1];
        double step = (x[n - 1] - x[0]) / m;
        for (int i = 0; i <= m; i++){
            pts[i] = x[0] + step * i;
        }

        var outCsv = new List<string> { "x,InterpolatedY" };
        foreach (var z in pts)
        {
            double v = spline.Evaluate(z);
            Console.WriteLine($"{z:F4}\t{v:F4}");
            outCsv.Add(
                z.ToString(CultureInfo.InvariantCulture)
                + ","
                + v.ToString(CultureInfo.InvariantCulture)
            );
        }

        File.WriteAllLines("evaluations.csv", outCsv);
    }
}
