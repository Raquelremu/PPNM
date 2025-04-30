using System;
using System.Collections;
using System.Collections.Generic;
using static System.Math;

public class ode {
    // Vector class implementation
    public class vector {
        private double[] data;
        
        public int size => data.Length;
        public double this[int i] {
            get => data[i];
            set => data[i] = value;
        }
        
        public vector(int n) { data = new double[n]; }
        public vector(params double[] values) { data = (double[])values.Clone(); }
        
        public static vector operator+(vector a, vector b) {
            vector result = new vector(a.size);
            for(int i=0; i<a.size; i++) result[i] = a[i] + b[i];
            return result;
        }
        
        public static vector operator-(vector a, vector b) {
            vector result = new vector(a.size);
            for(int i=0; i<a.size; i++) result[i] = a[i] - b[i];
            return result;
        }
        
        public static vector operator*(vector a, double x) {
            vector result = new vector(a.size);
            for(int i=0; i<a.size; i++) result[i] = a[i] * x;
            return result;
        }
        
        public static vector operator*(double x, vector a) => a * x;
        
        public static vector operator/(vector a, double x) {
            vector result = new vector(a.size);
            for(int i=0; i<a.size; i++) result[i] = a[i] / x;
            return result;
        }
        
        public double norm() {
            double sum = 0;
            foreach(double val in data) sum += val * val;
            return Sqrt(sum);
        }
        
        public vector copy() {
            return new vector((double[])data.Clone());
        }
    }

    // Generic list implementation
    public class genlist<T> : IEnumerable<T> {
        private T[] data;
        public int size => data.Length;
        public T this[int i] {
            get => data[i];
            set => data[i] = value;
        }
        
        public genlist() { data = Array.Empty<T>(); }
        
        public void add(T item) {
            T[] newdata = new T[size + 1];
            Array.Copy(data, newdata, size);
            newdata[size] = item;
            data = newdata;
        }
        
        public IEnumerator<T> GetEnumerator() {
            for(int i=0; i<size; i++) yield return data[i];
        }
        
        IEnumerator IEnumerable.GetEnumerator() {
            return GetEnumerator();
        }
    }

    // Embedded Runge-Kutta stepper (RK12: midpoint-Euler method)
    public static (vector, vector) rkstep12(
    	Func<double, vector, vector> f,
	double x,
	vector y,
	double h
   ) {
	    vector k0 = f(x, y);
	    vector k1 = f(x + h/2, y + k0*(h/2));
	    vector yh = y + k1*h;
	    vector δy = (k1 - k0)*h;
	    
	    // Add NaN/infinity check
	    for (int i = 0; i < yh.size; i++) {
		if (double.IsNaN(yh[i]) || double.IsInfinity(yh[i])) {
		    throw new Exception("Numerical instability in rkstep12");
		}
	    }
    
    return (yh, δy);
    }

    // Embedded Runge-Kutta stepper (RK23: Bogacki-Shampine method)
    public static (vector, vector) rkstep23(
    Func<double, vector, vector> f,
    double x,
    vector y,
    double h
) {
    vector k0 = f(x, y);
    vector k1 = f(x + h/2, y + k0*(h/2));
    vector k2 = f(x + 3*h/4, y + k1*(3*h/4));
    vector k3 = f(x + h, y + k0*(2*h/9) + k1*(h/3) + k2*(4*h/9));
    
    vector yh = y + k0*(2*h/9) + k1*(h/3) + k2*(4*h/9);
    vector yh_lower = y + k0*(7*h/24) + k1*(h/4) + k2*(h/3) + k3*(h/8);
    vector δy = yh - yh_lower;
    
    // Add NaN/infinity check
    for (int i = 0; i < yh.size; i++) {
        if (double.IsNaN(yh[i]) || double.IsInfinity(yh[i])) {
            throw new Exception("Numerical instability in rkstep23");
        }
    }
    
    return (yh, δy);
}

    // Adaptive-step-size driver routine
public static (genlist<double>, genlist<vector>) driver(
    Func<double, vector, vector> f,  // the f from dy/dx=f(x,y)
    (double, double) interval,      // (initial-point, final-point)
    vector yinit,                   // y(initial-point)
    double h = 0.125,               // initial step-size
    double acc = 0.01,              // absolute accuracy goal
    double eps = 0.01,              // relative accuracy goal
    string method = "rk12"          // integration method
) {
    var (a, b) = interval;
    double x = a;
    vector y = yinit.copy();
    var xlist = new genlist<double>(); xlist.add(x);
    var ylist = new genlist<vector>(); ylist.add(y.copy());
    
    // Add step counter for debug
    int stepCount = 0;
    const int maxSteps = 1000000;
    
    do {
        stepCount++;
        if (stepCount % 1000 == 0) {
            Console.WriteLine($"Step {stepCount}: x = {x}, h = {h}");
        }
        
        if (stepCount > maxSteps) {
            Console.WriteLine("Warning: Maximum step count reached");
            return (xlist, ylist);
        }
        
        if (x >= b) return (xlist, ylist); // job done
        if (x + h > b) h = b - x;         // last step should end at b
        
        (vector yh, vector δy) = method == "rk23" 
            ? rkstep23(f, x, y, h) 
            : rkstep12(f, x, y, h);
        
        // Add NaN/infinity check
        for (int i = 0; i < yh.size; i++) {
            if (double.IsNaN(yh[i]) || double.IsInfinity(yh[i])) {
                Console.WriteLine($"Numerical instability detected at x = {x}, h = {h}");
                return (xlist, ylist);
            }
        }
        
        double tol = (acc + eps * yh.norm()) * Sqrt(h/(b-a));
        double err = δy.norm();
        
        if (err <= tol) { // accept step
            x += h; y = yh;
            xlist.add(x);
            ylist.add(y.copy());
        }
        
        // Adjust step-size
        if (err > 0) h *= Min(Pow(tol/err, 0.25) * 0.95, 2);
        else h *= 2;
    } while (true);
}

    // Quadratic spline interpolation
    public static Func<double, vector> make_qspline(genlist<double> x, genlist<vector> y) {
        int n = x.size;
        vector[] b = new vector[n-1];
        vector[] c = new vector[n-1];
        
        // Calculate b_i and c_i coefficients
        for (int i = 0; i < n-1; i++) {
            double h = x[i+1] - x[i];
            vector dy = y[i+1] - y[i];
            b[i] = dy/h;
            if (i > 0) {
                c[i] = (b[i] - b[i-1])/(2*h);
            }
        }
        c[0] = new vector(y[0].size); // Initialize first c to zero vector
        
        return delegate(double z) {
            int i = binsearch(x, z);
            double dx = z - x[i];
            return y[i] + dx * (b[i] + c[i] * dx);
        };
    }

    // Binary search for spline
    private static int binsearch(genlist<double> x, double z) {
        int i = 0, j = x.size - 1;
        while (j - i > 1) {
            int mid = (i + j)/2;
            if (z > x[mid]) i = mid;
            else j = mid;
        }
        return i;
    }

    // Alternative interface with quadratic spline
    public static Func<double, vector> make_ode_ivp_qspline(
        Func<double, vector, vector> f,
        (double, double) interval,
        vector yinit,
        double acc = 0.01,
        double eps = 0.01,
        double hstart = 0.01,
        string method = "rk12"
    ) {
        var (xlist, ylist) = driver(f, interval, yinit, hstart, acc, eps, method);
        return make_qspline(xlist, ylist);
    }
}
