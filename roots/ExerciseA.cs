using System;
using static System.Math;
using System.Collections.Generic;

public class vector : List<double> {
    // Constructors
    public vector() : base() {}
    public vector(int n) : base(n) { for(int i=0; i<n; i++) Add(0.0); }
    public vector(vector v) : base(v) {}
    public vector(IEnumerable<double> d) : base(d) {}
    
    // Vector operations
    public static vector operator+(vector v, vector u) {
        vector r = new vector(v.Count);
        for(int i=0; i<v.Count; i++) r[i] = v[i] + u[i];
        return r;
    }
    
    // Unary minus operator (FIX for the error)
    public static vector operator-(vector v) {
        vector r = new vector(v.Count);
        for(int i=0; i<v.Count; i++) r[i] = -v[i];
        return r;
    }
    
    public static vector operator-(vector v, vector u) {
        vector r = new vector(v.Count);
        for(int i=0; i<v.Count; i++) r[i] = v[i] - u[i];
        return r;
    }
    
    public static vector operator*(vector v, double a) {
        vector r = new vector(v.Count);
        for(int i=0; i<v.Count; i++) r[i] = v[i] * a;
        return r;
    }
    
    public static vector operator*(double a, vector v) { return v*a; }
    
    public double norm() {
        double sum = 0;
        foreach(double vi in this) sum += vi*vi;
        return Sqrt(sum);
    }
    
    public vector copy() { return new vector(this); }
    
    public void map(Func<double, double> f) {
        for(int i=0; i<Count; i++) this[i] = f(this[i]);
    }
}

public class matrix {
    private int n, m;
    private double[,] data;
    
    public matrix(int n, int m) {
        this.n = n;
        this.m = m;
        data = new double[n,m];
    }
    
    public matrix(vector[] cols) {
        n = cols[0].Count;
        m = cols.Length;
        data = new double[n,m];
        for(int i=0; i<n; i++)
            for(int j=0; j<m; j++)
                data[i,j] = cols[j][i];
    }
    
    public double this[int i, int j] {
        get { return data[i,j]; }
        set { data[i,j] = value; }
    }
    
    public int size1 { get { return n; } }
    public int size2 { get { return m; } }
    
    public vector solve(vector b) {
        if(n != m) throw new ArgumentException("Matrix must be square");
        if(n != b.Count) throw new ArgumentException("Matrix and vector dimensions must match");
        
        // Gaussian elimination with partial pivoting
        matrix A = new matrix(n, n);
        vector x = new vector(n);
        for(int i=0; i<n; i++) {
            x[i] = b[i];
            for(int j=0; j<n; j++) A[i,j] = data[i,j];
        }
        
        for(int i=0; i<n-1; i++) {
            // Partial pivoting
            int maxrow = i;
            for(int k=i+1; k<n; k++)
                if(Abs(A[k,i]) > Abs(A[maxrow,i])) maxrow = k;
            
            if(maxrow != i) {
                for(int j=i; j<n; j++) {
                    double temp = A[i,j];
                    A[i,j] = A[maxrow,j];
                    A[maxrow,j] = temp;
                }
                double temp2 = x[i];
                x[i] = x[maxrow];
                x[maxrow] = temp2;
            }
            
            // Elimination
            for(int k=i+1; k<n; k++) {
                double factor = A[k,i]/A[i,i];
                for(int j=i; j<n; j++) A[k,j] -= factor * A[i,j];
                x[k] -= factor * x[i];
            }
        }
        
        // Back substitution
        for(int i=n-1; i>=0; i--) {
            double sum = 0;
            for(int j=i+1; j<n; j++) sum += A[i,j] * x[j];
            x[i] = (x[i] - sum) / A[i,i];
        }
        
        return x;
    }
}

public static class RootFinder {
    public static matrix jacobian(Func<vector, vector> f, vector x, vector fx=null, vector dx=null) {
        if(dx == null) {
            dx = x.copy();
            dx.map(xi => Max(Abs(xi), 1) * Pow(2, -26));
        }
        if(fx == null) fx = f(x);
        
        matrix J = new matrix(x.Count, x.Count);
        for(int j=0; j<x.Count; j++) {
            x[j] += dx[j];
            vector df = f(x) - fx;
            for(int i=0; i<x.Count; i++) J[i,j] = df[i]/dx[j];
            x[j] -= dx[j];
        }
        return J;
    }
    
    public static vector newton(
        Func<vector, vector> f,
        vector start,
        double acc=1e-2,
        vector δx=null
    ) {
        vector x = start.copy();
        vector fx = f(x), z, fz;
        double λmin = 1.0/128;
        
        do {
            if(fx.norm() < acc) break;
            
            matrix J = jacobian(f, x, fx, δx);
            vector Dx = J.solve(-fx); // Uses the unary minus operator
            
            double λ = 1;
            do {
                z = x + λ * Dx;
                fz = f(z);
                if(fz.norm() < (1-λ/2)*fx.norm()) break;
                if(λ < λmin) break;
                λ /= 2;
            } while(true);
            
            x = z; fx = fz;
        } while(true);
        
        return x;
    }
}

class ExerciseA {
    static vector rosenbrock_gradient(vector v) {
        double x = v[0], y = v[1];
        return new vector {
            2*(x-1) + 400*x*(x*x - y),
            200*(y - x*x)
        };
    }
    
    static vector himmelblau_gradient(vector v) {
        double x = v[0], y = v[1];
        return new vector {
            4*x*(x*x + y - 11) + 2*(x + y*y - 7),
            2*(x*x + y - 11) + 4*y*(x + y*y - 7)
        };
    }
    
    static void Main() {
    
        vector rosenbrock_root = RootFinder.newton(rosenbrock_gradient, new vector {0.5, 0.5});
        Console.WriteLine("Rosenbrock extremum at: ({0}, {1})", rosenbrock_root[0], rosenbrock_root[1]);
        
        // Test Newton's method on Himmelblau's function
        vector himmelblau_root1 = RootFinder.newton(himmelblau_gradient, new vector {3.0, 2.0});
        Console.WriteLine("Himmelblau minimum 1 at: ({0}, {1})", himmelblau_root1[0], himmelblau_root1[1]);
        
        vector himmelblau_root2 = RootFinder.newton(himmelblau_gradient, new vector {-2.0, 3.0});
        Console.WriteLine("Himmelblau minimum 2 at: ({0}, {1})", himmelblau_root2[0], himmelblau_root2[1]);
        
        vector himmelblau_root3 = RootFinder.newton(himmelblau_gradient, new vector {-3.0, -3.0});
        Console.WriteLine("Himmelblau minimum 3 at: ({0}, {1})", himmelblau_root3[0], himmelblau_root3[1]);
        
        vector himmelblau_root4 = RootFinder.newton(himmelblau_gradient, new vector {3.0, -2.0});
        Console.WriteLine("Himmelblau minimum 4 at: ({0}, {1})", himmelblau_root4[0], himmelblau_root4[1]);
    }
}
