using System;

class Program
{
    static void Main()
    {
        for (int i = 1; i <= 10; i++)
        {
            Console.WriteLine($"Γ({i}) = {sfuns.fgamma(i)}");
        }
    }
}
