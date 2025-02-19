using System;

public class vec
{
    public double x, y, z;

    // Default constructor
    public vec() { x = y = z = 0; }

    // Parameterized constructor
    public vec(double x, double y, double z) { this.x = x; this.y = y; this.z = z; }

    // Override ToString method for easy printing
    public override string ToString() => $"{x} {y} {z}";

    // Print method for debugging
    public void print(string s = "") { Console.WriteLine($"{s}{x} {y} {z}"); }

    // Operator Overloads
    public static vec operator *(vec v, double c) => new vec(c * v.x, c * v.y, c * v.z);
    public static vec operator *(double c, vec v) => v * c;
    public static vec operator /(vec v, double c) => new vec(v.x / c, v.y / c, v.z / c);
    public static vec operator +(vec u, vec v) => new vec(u.x + v.x, u.y + v.y, u.z + v.z);
    public static vec operator -(vec u) => new vec(-u.x, -u.y, -u.z);
    public static vec operator -(vec u, vec v) => new vec(u.x - v.x, u.y - v.y, u.z - v.z);

    // Dot Product
    public double dot(vec other) => this.x * other.x + this.y * other.y + this.z * other.z;
    public static double dot(vec v, vec w) => v.x * w.x + v.y * w.y + v.z * w.z;

    // Cross Product
    public vec cross(vec v) =>
        new vec(
            y * v.z - z * v.y,
            z * v.x - x * v.z,
            x * v.y - y * v.x
        );

    // Norm
    public double norm() => Math.Sqrt(x * x + y * y + z * z);

    // Approximate comparison
    static bool approx(double a, double b, double acc = 1e-9, double eps = 1e-9)
    {
        if (Math.Abs(a - b) < acc) return true;
        if (Math.Abs(a - b) < (Math.Abs(a) + Math.Abs(b)) * eps) return true;
        return false;
    }

    public bool approx(vec other)
    {
        if (!approx(this.x, other.x)) return false;
        if (!approx(this.y, other.y)) return false;
        if (!approx(this.z, other.z)) return false;
        return true;
    }

    public static bool approx(vec u, vec v) => u.approx(v);
}

