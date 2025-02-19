using System;

class Program
{
    static void Main()
    {
        vec v1 = new vec(1, 2, 3);
        vec v2 = new vec(4, 5, 6);

        Console.WriteLine("Vector v1: " + v1);
        Console.WriteLine("Vector v2: " + v2);

        // Addition
        Console.WriteLine("v1 + v2 = " + (v1 + v2));

        // Subtraction
        Console.WriteLine("v1 - v2 = " + (v1 - v2));

        // Scalar Multiplication
        Console.WriteLine("2 * v1 = " + (2 * v1));

        // Dot Product
        Console.WriteLine("Dot product v1 . v2 = " + v1.dot(v2));

        // Cross Product
        Console.WriteLine("Cross product v1 x v2 = " + v1.cross(v2));

        // Norm
        Console.WriteLine("Norm of v1 = " + v1.norm());

        // Approximate Comparison
        Console.WriteLine("Are v1 and v2 approximately equal? " + vec.approx(v1, v2));
    }
}

