using System;
using System.IO;
using System.Diagnostics;

public class vector {
    public double[] data;
    public int size => data.Length;
    public double this[int i] {
        get => data[i];
        set => data[i] = value;
    }
    public vector(int n) {
        data = new double[n];
    }
    public void Print(string name = "Vector") {
        Console.Write(name + ": [");
        for (int i = 0; i < size; i++) {
            Console.Write($"{data[i]:F12}");
            if (i < size - 1) Console.Write(", ");
        }
        Console.WriteLine("]");
    }
}

public class matrix {
    public readonly int size1, size2;
    private double[] data;
    public matrix(int n, int m) {
        size1 = n; size2 = m;
        data = new double[n * m];
    }
    public double M[int i, int j] {
        get => data[i + j * size1];
        set => data[i + j * size1] = value;
    }
    public void Print(string name = "Matrix") {
        Console.WriteLine(name + ":");
        for (int i = 0; i < size1; i++) {
            for (int j = 0; j < size2; j++) {
                Console.Write($"{this[i, j]:F12} ");
            }
            Console.WriteLine();
        }
    }
    public matrix Copy() {
        matrix copy = new matrix(size1, size2);
        for (int i = 0; i < size1; i++) {
            for (int j = 0; j < size2; j++) {
                copy[i, j] = M[i, j];
            }
        }
        return copy;
    }
    public static matrix Multiply(matrix A, matrix B) {
        if (A.size2 != B.size1) throw new Exception("Matrix dimensions do not match for multiplication");
        matrix C = new matrix(A.size1, B.size2);
        for (int i = 0; i < A.size1; i++) {
            for (int j = 0; j < B.size2; j++) {
                double sum = 0.0;
                for (int k = 0; k < A.size2; k++) {
                    sum += A[i, k] * B[k, j];
                }
                C[i, j] = sum;
            }
        }
        return C;
    }
    public matrix Transpose() {
        matrix trans = new matrix(size2, size1);
        for (int i = 0; i < size1; i++) {
            for (int j = 0; j < size2; j++) {
                trans[j, i] = this[i, j];
            }
        }
        return trans;
    }
}

public class QR {
    public matrix Q, R;
    public QR(matrix A) {
        int n = A.size1;
        int m = A.size2;
        Q = A.Copy();
        R = new matrix(m, m);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < i; j++) {
                double dot = 0.0;
                for (int k = 0; k < n; k++)
                    dot += Q[k, j] * Q[k, i];
                R[j, i] = dot;
                for (int k = 0; k < n; k++)
                    Q[k, i] -= dot * Q[k, j];
            }
            double norm = 0.0;
            for (int k = 0; k < n; k++)
                norm += Q[k, i] * Q[k, i];
            norm = Math.Sqrt(norm);
            R[i, i] = norm;
            for (int k = 0; k < n; k++)
                Q[k, i] /= norm;
        }
    }

    public vector solve(vector b, matrix Qt) {
        int n = Q.size1;
        int m = R.size1;
        if (b.size != n)
            throw new Exception("Vector b size does not match matrix dimensions");

        vector y = new vector(m);
        for (int i = 0; i < m; i++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++)
                sum += Qt[i, k] * b[k];
            y[i] = sum;
        }
        vector x = new vector(m);
        for (int i = m - 1; i >= 0; i--) {
            double sum = y[i];
            for (int k = i + 1; k < m; k++)
                sum -= R[i, k] * x[k];
            if (Math.Abs(R[i, i]) < 1e-14)
                throw new Exception($"Singular matrix at row {i}");
            x[i] = sum / R[i, i];
        }
        return x;
    }
    public matrix Inverse() {
        int n = Q.size1;
        matrix I = new matrix(n, n);
        matrix B = new matrix(n, n);
        for (int i = 0; i < n; i++)
            I[i, i] = 1.0;
        matrix Qt = new matrix(Q.size2, n);
        for (int i = 0; i < Q.size2; i++)
            for (int j = 0; j < n; j++)
                Qt[i, j] = Q[j, i];
        for (int i = 0; i < n; i++) {
            vector ei = new vector(n);
            for (int j = 0; j < n; j++)
                ei[j] = I[j, i];
            vector x = solve(ei, Qt);
            for (int j = 0; j < n; j++)
                B[j, i] = x[j];
        }
        return B;
    }
    public double det() {
        double determinant = 1.0;
        for (int i = 0; i < R.size1; i++)
            determinant *= R[i, i];
        return determinant;
    }
}

public class LeastSquares {
    public static (vector, matrix) lsfit(Func<double, double>[] fs, vector x, vector y, vector dy) {
        int n = x.size;
        int m = fs.Length;
        matrix A = new matrix(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                A[i, j] = fs[j](x[i]) / dy[i];
            }
        }
        vector b = new vector(n);
        for (int i = 0; i < n; i++) {
            b[i] = y[i] / dy[i];
        }
        
        QR qr = new QR(A);
        matrix Qt = new matrix(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                Qt[i, j] = qr.Q[j, i];
            }
        }
        vector coeff = qr.solve(b, Qt);
        
        matrix Rinv = qr.R.Inverse();
        matrix cov = matrix.Multiply(Rinv, Rinv.Transpose());
        
        return (coeff, cov);
    }
}

public class Program {
    public static void Main(string[] args) {
        double[] t_data = { 1, 2, 3, 4, 6, 9, 10, 13, 15 };
        double[] y_data = { 117, 100, 88, 72, 53, 29.5, 25.2, 15.2, 11.1 };
        double[] dy_data = { 6, 5, 4, 4, 4, 3, 3, 2, 2 };
        int n = t_data.Length;

        vector t = new vector(n);
        vector ln_y = new vector(n);
        vector dln_y = new vector(n);
        for (int i = 0; i < n; i++) {
            t[i] = t_data[i];
            ln_y[i] = Math.Log(y_data[i]);
            dln_y[i] = dy_data[i] / y_data[i];
        }

        Func<double, double>[] fs = {
            z => 1.0,
            z => -z
        };

        (vector coeff, matrix cov) = LeastSquares.lsfit(fs, t, ln_y, dln_y);
        double ln_a = coeff[0];
        double lambda = coeff[1];
        double a = Math.Exp(ln_a);
        double half_life = Math.Log(2) / lambda;
        double modern_half_life = 3.6319;
        
        double dln_a = Math.Sqrt(cov[0, 0]);
        double dlambda = Math.Sqrt(cov[1, 1]);
        double da = a * dln_a; 
        double dhalf_life = (Math.Log(2) / (lambda * lambda)) * dlambda; 
        
        Console.WriteLine("Fitted parameters:");
        Console.WriteLine($"a = {a:F4} ± {da:F4}");
        Console.WriteLine($"λ = {lambda:F4} ± {dlambda:F4}");
        Console.WriteLine($"Half-life T½ = {half_life:F4} ± {dhalf_life:F4} days");
        Console.WriteLine($"Modern value (224Ra) = {modern_half_life:F4} days");
        Console.WriteLine($"Difference = {Math.Abs(half_life - modern_half_life):F4} days ({Math.Abs(half_life - modern_half_life) / modern_half_life * 100:F2}%)");
        
        Console.WriteLine("\nCovariance matrix:");
        cov.Print();

        using (StreamWriter writer = new StreamWriter("plot_data.txt")) {
            writer.WriteLine("# Experimental data: t, y, dy");
            for (int i = 0; i < n; i++) {
                writer.WriteLine($"{t_data[i]} {y_data[i]} {dy_data[i]}");
            }
            writer.WriteLine();
            writer.WriteLine("# Fitted data: t, y_fit");
            double t_min = 0, t_max = t_data[n - 1] + 2;
            for (double ti = t_min; ti <= t_max; ti += 0.1) {
                double ln_y_fit = ln_a - lambda * ti;
                double y_fit = Math.Exp(ln_y_fit);
                writer.WriteLine($"{ti} {y_fit}");
            }
        }
        Console.WriteLine("\nData for plotted");
    }
}
