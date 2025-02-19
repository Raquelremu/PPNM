using static System.Console;
using System.Threading;

class main
{
    public class datum { public long start, stop; public double sum; }

    public static void harm(object obj)
    {
        datum d = (datum)obj;
        d.sum = 0;
        for (long i = d.start + 1; i <= d.stop; i++) d.sum += 1.0 / i;
    }

    public static int Main(string[] argv)
    {
        long argc = argv.Length;
        long nthreads = 1, nterms = (long)1e8;

        for (long i = 0; i < argc; i++)
        {
            string arg = argv[i];
            if (arg == "-threads" && i + 1 < argc) nthreads = long.Parse(argv[i + 1]);
            if (arg == "-terms" && i + 1 < argc) nterms = (long)double.Parse(argv[i + 1]);
        }

        if (nthreads <= 0)
        {
            Error.WriteLine("Error: Number of threads must be at least 1.");
            return 1;
        }

        var threads = new Thread[nthreads];
        var data = new datum[nthreads];

        for (long i = 0; i < nthreads; i++)
        {
            data[i] = new datum();  // 🔴 Fix: Initialize each element before use
            data[i].start = nterms / nthreads * i;
            data[i].stop = nterms / nthreads * (i + 1);
        }

        for (long i = 0; i < nthreads; i++)
        {
            threads[i] = new Thread(harm);
            threads[i].Start(data[i]);
        }

        foreach (var thread in threads) thread.Join();
        double sum = 0;
        foreach (var d in data) 
        {
            if (d != null) // 🔴 Fix: Ensure `d` is not null before accessing `.sum`
            {
                sum += d.sum;
            }
        }

        WriteLine($"sum={sum}");
        return 0;
    }
}

