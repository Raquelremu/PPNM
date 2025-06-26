using System;
using System.Collections.Generic;
using static System.Math;
using System.IO;

class Program
{
    static int LastNewtonSteps;

    static void Main(string[] args)
    {
        // Task A
        vector x0 = new vector(2); 
        x0[0] = 1; 
        x0[1] = 1;
        
        vector rosenbrockMin = Newton(Rosenbrock, x0);
        Console.WriteLine($"Rosenbrock minimum at: {rosenbrockMin[0]:F6}, {rosenbrockMin[1]:F6}");
        Console.WriteLine($"Steps taken: {LastNewtonSteps}");
        Console.WriteLine($"Function value at minimum: {Rosenbrock(rosenbrockMin):F6}");
	Console.WriteLine($" ");
	
        x0[0] = -3; 
        x0[1] = 3;
        vector himmelblauMin = Newton(Himmelblau, x0);
        Console.WriteLine($"HimmelBlau minimum at: {himmelblauMin[0]:F6}, {himmelblauMin[1]:F6}");
        Console.WriteLine($"Steps taken: {LastNewtonSteps}");
        Console.WriteLine($"Function value at minimum: {Himmelblau(himmelblauMin):F6}");
	Console.WriteLine($" ");
	
        // Task B
        var energy = new List<double>();
        var signal = new List<double>();
        var error = new List<double>();

        string line;
        var separators = new char[] { ' ', '\t' };
        var options = StringSplitOptions.RemoveEmptyEntries;
        while ((line = Console.ReadLine()) != null)
        {
            line = line.Trim();
            if (line.Length == 0 || line.StartsWith("#")) continue;
            string[] parts = line.Split(separators, options);
            energy.Add(double.Parse(parts[0]));
            signal.Add(double.Parse(parts[1]));
            error.Add(double.Parse(parts[2]));
        }

        vector initialParams = new vector(3);
        initialParams[0] = 125.3;    
        initialParams[1] = 0.00407;  
        initialParams[2] = 6.0;      
        vector fitParams = Newton(p => Deviation(p, energy, signal, error), initialParams);
        Console.Write("Higgs fit parameters: "); fitParams.Print();

        using (var w1 = new StreamWriter("exp_data.txt"))
        {
            for (int i = 0; i < energy.Count; i++)
                w1.WriteLine($"{energy[i]} {signal[i]} {error[i]}");
        }
        using (var w2 = new StreamWriter("fit_curve.txt"))
        {
            double minE = energy[0], maxE = energy[energy.Count - 1];
            for (double E = minE; E <= maxE; E += 0.5)
                w2.WriteLine($"{E} {BreitWigner(E, fitParams[0], fitParams[1], fitParams[2])}");
        }
    }

    static vector Newton(Func<vector, double> φ, vector x, double acc = 1e-3, int maxSteps = 1000)
    {
        LastNewtonSteps = 0;
        for (int k = 0; k < maxSteps; k++)
        {
            LastNewtonSteps++;
            vector g = Gradient(φ, x);
            if (g.norm() < acc) break;
            matrix H = Hessian(φ, x);
            QR qrH = new QR(H);
            matrix Qt = new matrix(H.size2, H.size1);
            for (int i = 0; i < H.size2; i++)
                for (int j = 0; j < H.size1; j++)
                    Qt[i, j] = qrH.Q[j, i];
            vector dx = qrH.solve(-g, Qt);

            double λ = 1;
            while (λ > 1.0 / 128)
            {
                vector x_new = x.copy();
                for (int i = 0; i < x.size; i++) x_new[i] += λ * dx[i];
                if (φ(x_new) < φ(x)) break;
                λ /= 2;
            }
            for (int i = 0; i < x.size; i++) x[i] += λ * dx[i];
        }
        return x;
    }

    static vector Gradient(Func<vector, double> φ, vector x)
    {
    double φx = φ(x);
    vector g = new vector(x.size);
    for (int i = 0; i < x.size; i++)
    {
        double dxi = Abs(x[i]) * Pow(2, -26);
        x[i] += dxi;
        g[i] = (φ(x) - φx)/dxi;
        x[i] -= dxi;
    }
    return g;
    }

    static matrix Hessian(Func<vector, double> φ, vector x)
    {
    int n = x.size;
    matrix H = new matrix(n,n);
    vector g0 = Gradient(φ, x);
    double eps2 = Pow(2, -17);
    for (int j = 0; j < n; j++)
    {
        double dxj = Abs(x[j]) * eps2;
        x[j] += dxj;
        vector g1 = Gradient(φ, x);
        for (int i = 0; i < n; i++)
            H[i,j] = (g1[i] - g0[i]) / dxj;
        x[j] -= dxj;
    }
    return H;
    }

    static double Rosenbrock(vector v)
    {
        double x = v[0], y = v[1];
        return Pow(1 - x, 2) + 100 * Pow(y - x * x, 2);
    }

    static double Himmelblau(vector v)
    {
        double x = v[0], y = v[1];
        return Pow(x * x + y - 11, 2) + Pow(x + y * y - 7, 2);
    }

    static double BreitWigner(double E, double m, double Γ, double A)
        => A / (Pow(E - m, 2) + Γ * Γ / 4);

    static double Deviation(vector p, List<double> E, List<double> o, List<double> Δ)
    {
        double sum = 0;
        for (int i = 0; i < E.Count; i++)
        {
            double diff = (BreitWigner(E[i], p[0], p[1], p[2]) - o[i]) / Δ[i];
            sum += diff * diff;
        }
        return sum;
    }
}

public class vector
{
    public double[] data;
    public int size => data.Length;
    public double this[int i] { get => data[i]; set => data[i] = value; }
    public vector(int n) => data = new double[n];
    public vector(double[] d) => data = (double[])d.Clone();

    public double norm()
    {
        double sum = 0;
        foreach (var x in data) sum += x * x;
        return Sqrt(sum);
    }

    public vector copy()
    {
        var v = new vector(size);
        Array.Copy(data, v.data, size);
        return v;
    }

    public void Print(string s = "")
    {
        Console.Write(s);
        foreach (var x in data) Console.Write($"{x:F6} ");
        Console.WriteLine();
    }

    public static vector operator +(vector a, vector b)
    {
        var r = new vector(a.size);
        for (int i = 0; i < a.size; i++) r[i] = a[i] + b[i];
        return r;
    }

    public static vector operator -(vector a, vector b)
    {
        var r = new vector(a.size);
        for (int i = 0; i < a.size; i++) r[i] = a[i] - b[i];
        return r;
    }

    public static vector operator -(vector a)
    {
        var r = new vector(a.size);
        for (int i = 0; i < a.size; i++) r[i] = -a[i];
        return r;
    }

    public static vector operator*(vector a, double s)
    {
        var r = new vector(a.size);
        for (int i = 0; i < a.size; i++) r[i] = a[i] * s;
        return r;
    }

    public static vector operator*(double s, vector a) => a * s;
}

public class matrix
{
    public int size1, size2;
    private double[] data;
    public matrix(int n, int m)
    {
        size1 = n; size2 = m;
        data = new double[n * m];
    }
    public double this[int i, int j] { get => data[i + j * size1]; set => data[i + j * size1] = value; }
    public matrix copy()
    {
        var M = new matrix(size1, size2);
        Array.Copy(data, M.data, data.Length);
        return M;
    }
    public static matrix operator *(matrix A, matrix B)
    {
        if (A.size2 != B.size1) throw new Exception("Dimension mismatch");
        var C = new matrix(A.size1, B.size2);
        for (int i = 0; i < A.size1; i++)
            for (int j = 0; j < B.size2; j++)
                for (int k = 0; k < A.size2; k++)
                    C[i, j] += A[i, k] * B[k, j];
        return C;
    }
}

public class QR
{
    public matrix Q, R;
    public QR(matrix A)
    {
        int n = A.size1, m = A.size2;
        Q = A.copy(); R = new matrix(m, m);
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < i; j++)
            {
                double dot = 0;
                for (int k = 0; k < n; k++) dot += Q[k, j] * Q[k, i];
                R[j, i] = dot;
                for (int k = 0; k < n; k++) Q[k, i] -= dot * Q[k, j];
            }
            double norm = 0;
            for (int k = 0; k < n; k++) norm += Q[k, i] * Q[k, i];
            R[i, i] = Sqrt(norm);
            for (int k = 0; k < n; k++) Q[k, i] /= R[i, i];
        }
    }
    public vector solve(vector b, matrix Qt)
    {
        int n = Q.size1, m = R.size1;
        var y = new vector(m);
        for (int i = 0; i < m; i++)
            for (int k = 0; k < n; k++)
                y[i] += Qt[i, k] * b[k];
        var x = new vector(m);
        for (int i = m - 1; i >= 0; i--)
        {
            double sum = y[i];
            for (int k = i + 1; k < m; k++) sum -= R[i, k] * x[k];
            x[i] = sum / R[i, i];
        }
        return x;
    }
}

