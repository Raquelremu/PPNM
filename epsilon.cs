using System;

class Program
{
    static void Main()
    {
        // 1. Maximum and Minimum Representable Integers
        Console.WriteLine("### Maximum and Minimum Representable Integers ###");
        
        int maxInt = 1;
        while (maxInt + 1 > maxInt) 
            maxInt++;

        int minInt = -1;
        while (minInt - 1 < minInt) 
            minInt--;

        Console.WriteLine($"Maximum int: {maxInt}, Expected: {int.MaxValue}");
        Console.WriteLine($"Minimum int: {minInt}, Expected: {int.MinValue}");

// 2. Machine Epsilon Calculation
        Console.WriteLine("\n### Machine Epsilon Calculation ###");

        double epsilon = 1.0;
        while (1.0 + epsilon != 1.0)
            epsilon /= 2;
        epsilon *= 2; // The last valid epsilon before breaking the condition

        float epsilonFloat = 1f;
        while (1f + epsilonFloat != 1f)
            epsilonFloat /= 2;
        epsilonFloat *= 2; // The last valid epsilon before breaking the condition

        Console.WriteLine($"Double Machine Epsilon: {epsilon}, Expected: {Math.Pow(2, -52)}");
        Console.WriteLine($"Float Machine Epsilon: {epsilonFloat}, Expected: {Math.Pow(2, -23)}");

// 3. Tiny value comparison
        Console.WriteLine("\n### Tiny Value Comparison ###");

        double tiny = epsilon / 2;
        double a = 1.0 + tiny;
        double b = tiny + tiny + 1.0;

        Console.WriteLine($"a == b: {a == b}");
        Console.WriteLine($"a > 1: {a > 1}");
        Console.WriteLine($"b > 1: {b > 1}");

// 4. Comparing Doubles Using an Approximation Function
        Console.WriteLine("\n### Comparing Doubles with Approximation ###");

        double d1 = 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1;
        double d2 = 8 * 0.1;

        Console.WriteLine($"d1: {d1:e15}, d2: {d2:e15}");
        Console.WriteLine($"d1 == d2: {d1 == d2}");
        Console.WriteLine($"approx(d1, d2): {Approx(d1, d2)}");
    }

    static bool Approx(double a, double b, double acc = 1e-9, double eps = 1e-9)
    {
        if (Math.Abs(b - a) <= acc) return true;
        if (Math.Abs(b - a) <= Math.Max(Math.Abs(a), Math.Abs(b)) * eps) return true;
        return false;
    }
}

