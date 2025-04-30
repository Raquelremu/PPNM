using System;
using static System.Math;

public class Program {

    public class vector {
        public double[] data;
        public int size { get { return data.Length; } }

        public vector(int n) {
            data = new double[n];
        }

        public double this[int i] {
            get { return data[i]; }
            set { data[i] = value; }
        }
    }

    public class matrix {
        public double[][] data;
        public int size1 { get { return data.Length; } }
        public int size2 { get { return data[0].Length; } }

        public matrix(int n, int m) {
            data = new double[n][];
            for (int i = 0; i < n; i++){
                data[i] = new double[m];
            }
        }

        public matrix(int n) : M(n, n) {}

        public double this[int i, int j] {
            get { return data[i][j]; }
            set { data[i][j] = value; }
        }

        public matrix copy(){
            matrix Copy = new matrix(size1, size2);
            for (int i = 0; i < size1; i++){
                for (int j = 0; j < size2; j++){
                    Copy[i, j] = this[i, j];
                }
            }
            return Copy;
        }

        public static matrix id(int n){
            matrix I = new matrix(n, n);
            for (int i = 0; i < n; i++){
                I[i, i] = 1;
            }
            return I;
        }
    }

    public static void timesJ(matrix A, int p, int q, double theta){
        double c = Cos(theta), s = Sin(theta);
        for (int i = 0; i < A.size1; i++){
            double aip = A[i, p];
            double aiq = A[i, q];
            A[i, p] = c * aip - s * aiq;
            A[i, q] = s * aip + c * aiq;
        }
    }

    public static void Jtimes(matrix A, int p, int q, double theta){
        double c = Cos(theta), s = Sin(theta);
        for (int j = 0; j < A.size2; j++){
            double apj = A[p, j];
            double aqj = A[q, j];
            A[p, j] =  c * apj + s * aqj;
            A[q, j] = -s * apj + c * aqj;
        }
    }

    public static (vector, matrix) cyclic(matrix M){
        matrix A = M.copy();
        int n = A.size1;
        matrix V = matrix.id(n);
        vector w = new vector(n); 

        bool changed;
        do {
            changed = false;
            for (int p = 0; p < n - 1; p++){
                for (int q = p + 1; q < n; q++){
                    double app = A[p, p];
                    double aqq = A[q, q];
                    double apq = A[p, q];
                    double theta = 0.5 * Atan2(2 * apq, aqq - app);
                    double c = Cos(theta), s = Sin(theta);

                    double new_app = c * c * app - 2 * s * c * apq + s * s * aqq;
                    double new_aqq = s * s * app + 2 * s * c * apq + c * c * aqq;

                    if (new_app != app || new_aqq != aqq) {
                        changed = true;
 
                        timesJ(A, p, q, theta);

                        Jtimes(A, p, q, -theta);

                        timesJ(V, p, q, theta);
                    }
                }
            }
        } while(changed);

        for (int i = 0; i < n; i++){
            w[i] = A[i, i];
        }

        return (w, V);
    }

    public static matrix Transpose(matrix A) {
        matrix T = new matrix(A.size2, A.size1);
        for (int i = 0; i < A.size1; i++){
            for (int j = 0; j < A.size2; j++){
                T[j, i] = A[i, j];
            }
        }
        return T;
    }

    public static matrix Multiply(matrix A, matrix B) {
        if (A.size2 != B.size1) {
            throw new Exception("Multiplication not possible");
        }
        matrix C = new matrix(A.size1, B.size2);
        for (int i = 0; i < A.size1; i++){
            for (int j = 0; j < B.size2; j++){
                double sum = 0;
                for (int k = 0; k < A.size2; k++){
                    sum += A[i, k] * B[k, j];
                }
                C[i, j] = sum;
            }
        }
        return C;
    }

    public static matrix Diag(vector v) {
        matrix D = new matrix(v.size, v.size);
        for (int i = 0; i < v.size; i++){
            D[i, i] = v[i];
        }
        return D;
    }

    public static void PrintMatrix(matrix A) {
        for (int i = 0; i < A.size1; i++){
            for (int j = 0; j < A.size2; j++){
                Console.Write($"{A[i, j],10:F5} ");
            }
            Console.WriteLine();
        }
        Console.WriteLine();
    }

    static void Main(string[] args){

        double rmax = 10;
        double dr   = 0.3;

        for (int i = 0; i < args.Length; i++){
            if (args[i] == "-rmax") rmax = double.Parse(args[i+1]);
            if (args[i] == "-dr")   dr   = double.Parse(args[i+1]);
        }

        int npoints = (int)(rmax / dr) - 1;
        vector r = new vector(npoints);
        for (int i = 0; i < npoints; i++){
            r[i] = dr * (i + 1);
        }

        matrix H = new matrix(npoints, npoints);
        double diagval    = -2.0 * (-0.5 / (dr * dr));
        double offdiagval =  1.0 * (-0.5 / (dr * dr));   
        for (int i = 0; i < npoints - 1; i++){
            H[i, i]   = diagval;
            H[i, i+1] = offdiagval;
            H[i+1, i] = offdiagval;
        }
        H[npoints - 1, npoints - 1] = diagval;

        for (int i = 0; i < npoints; i++){
            H[i, i] += -1.0 / r[i];
        }

        var (e, V) = cyclic(H);

        Console.WriteLine($"# rmax={rmax}, dr={dr}, npoints={npoints}");
        Console.WriteLine("# Lowest eigenvalues (s-wave energies)");
        for (int k = 0; k < 5 && k < e.size; k++){
            Console.WriteLine($"E[{k}] = {e[k]}");
        }
        Console.WriteLine();

        matrix D = Diag(e);

        matrix Vt = Transpose(V);
        matrix VtHV = Multiply(Multiply(Vt, H), V);
        Console.WriteLine("Matrix V^T * H * V");
        PrintMatrix(VtHV);

        matrix VDVt = Multiply(Multiply(V, D), Transpose(V));
        Console.WriteLine("Matrix V * D * V^T");
        PrintMatrix(VDVt);

        matrix VVt = Multiply(V, Transpose(V));
        Console.WriteLine("Matrix V * V^T ");
        PrintMatrix(VVt);

        matrix VtV = Multiply(Transpose(V), V);
        Console.WriteLine("Matrix V^T * V");
        PrintMatrix(VtV);
    }
}
