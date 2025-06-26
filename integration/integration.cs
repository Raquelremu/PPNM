using System;
using System.Diagnostics;
using System.IO;
using static System.Math;

public class AdaptiveIntegration
{
    public static int FunctionEvaluations = 0;

    public static (double result, double error) integrate(Func<double, double> f, double a, double b,
        double acc = 0.001, double eps = 0.001, double f2 = double.NaN, double f3 = double.NaN)
    {
        double h = b - a;
        if (double.IsNaN(f2)) { f2 = CountedEval(f, a + 2 * h / 6); f3 = CountedEval(f, a + 4 * h / 6); }
        double f1 = CountedEval(f, a + h / 6), f4 = CountedEval(f, a + 5 * h / 6);
        double Q = (2 * f1 + f2 + f3 + 2 * f4) / 6 * (b - a);
        double q = (f1 + f2 + f3 + f4) / 4 * (b - a);
        double err = Abs(Q - q);
        if (err <= acc + eps * Abs(Q)) return (Q, err);
        else 
        {
            var (Q1, err1) = integrate(f, a, (a + b) / 2, acc / Sqrt(2), eps, f1, f2);
            var (Q2, err2) = integrate(f, (a + b) / 2, b, acc / Sqrt(2), eps, f3, f4);
            return (Q1 + Q2, Sqrt(err1 * err1 + err2 * err2));
        }
    }

    public static double Erf(double z)
    {
        if (z < 0) return -Erf(-z);
        if (z <= 1) return 2 / Sqrt(PI) * integrate(x => Exp(-x * x), 0, z).result;
        return 1 - 2 / Sqrt(PI) * integrate(t => Exp(-Pow(z + (1 - t) / t, 2)) / (t * t), 0, 1).result;
    }

    public static void GenerateErrorAnalysisData()
    {
        using (var writer = new StreamWriter("error_data.dat"))
        {
            writer.WriteLine("acc estimated_err actual_err evaluations");
            double exact = 0.84270079294971486934;
            for (int i = 1; i <= 10; i++)
            {
                double acc = Pow(10, -i);
                FunctionEvaluations = 0;
                var (result, estErr) = integrate(x => 2/Sqrt(PI)*Exp(-x*x), 0, 1, acc, 0);
                writer.WriteLine($"{acc:E5} {estErr:E5} {Abs(result - exact):E5} {FunctionEvaluations}");
            }
        }
    }

    public static void GenerateComparisonData()
    {
        var testCases = new (string name, Func<double, double> f, double a, double b, double exact)[]
        {
            ("1/sqrt(x)", x => 1/Sqrt(x), 0, 1, 2.0),
            ("ln(x)/sqrt(x)", x => Log(x)/Sqrt(x), 0, 1, -4.0),
            ("exp(-x^2)", x => Exp(-x*x), 0, 1, Sqrt(PI)/2)
        };

        using (var writer = new StreamWriter("comparison_results.dat"))
        {
            writer.WriteLine("# Function C#Result ExactValue C#Error");
            foreach (var test in testCases)
            {
                FunctionEvaluations = 0;
                var (result, error) = integrate(test.f, test.a, test.b, 1e-6, 1e-6);
                
                writer.WriteLine($"{test.name} {result:F6} {test.exact:F6} {error:E3}");
            }
        }
    }

    public static void CreateGnuplotScripts()
    {
        File.WriteAllText("plot_errors.gp", 
            "set terminal pngcairo enhanced font 'Arial,12'\n" +
            "set output 'errors.png'\n" +
            "set logscale xy\n" +
            "set format x '10^{%L}'\n" +
            "set format y '10^{%L}'\n" +
            "set title 'Error Analysis (log-log scale)'\n" +
            "set xlabel 'Requested Accuracy'\n" +
            "set ylabel 'Error Magnitude'\n" +
            "set key top left\n" +
            "plot 'error_data.dat' using 1:2 with linespoints pt 7 title 'Estimated', " +
            "'error_data.dat' using 1:3 with linespoints pt 9 title 'Actual'");

        File.WriteAllText("plot_evaluations.gp", 
            "set terminal pngcairo enhanced font 'Arial,12'\n" +
            "set output 'evaluations.png'\n" +
            "set logscale xy\n" +
            "set format x '10^{%L}'\n" +
            "set format y '%g'\n" +
            "set title 'Function Evaluations vs Accuracy'\n" +
            "set xlabel 'Requested Accuracy'\n" +
            "set ylabel 'Evaluations'\n" +
            "plot 'error_data.dat' using 1:4 with linespoints pt 5 lw 2 title ''");

        File.WriteAllText("plot_comparison.gp", 
            "set terminal pngcairo enhanced font 'Arial,12'\n" +
            "set output 'comparison.png'\n" +
            "set style data histograms\n" +
            "set style fill solid\n" +
            "set boxwidth 0.8\n" +
            "set title 'Integration Results (with error estimates)'\n" +
            "set ylabel 'Value'\n" +
            "set xtics rotate by -45\n" +
            "set yrange [-5:3]\n" + 
            "plot 'comparison_results.dat' using 2:xtic(1) title 'C#', " +
            "'comparison_results.dat' using 3 title 'Exact', " +
            "'' using 0:($2>$3?$2:$3):4 with labels offset 0,1 font ',10' notitle");
    }

    private static double CountedEval(Func<double, double> f, double x)
    {
        FunctionEvaluations++;
        return f(x);
    }

    public static void Main(string[] args)
    {
        GenerateErrorAnalysisData();
        GenerateComparisonData(); 
        CreateGnuplotScripts();  

        Console.WriteLine("Generated all data files and plot scripts");
        Console.WriteLine("Run with: gnuplot plot_errors.gp, plot_evaluations.gp, plot_comparison.gp");
    }
}
