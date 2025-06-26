using System;
using System.Diagnostics;

public class vector{
   public double[] data;
   public int size => data.Length;
   public double this[int i]{
      get => data[i];
      set => data[i] = value;
   }
   public vector(int n){
      data = new double[n];
   }
   public void Print(string name = "Vector"){ 
      Console.Write(name + ": [");
      for (int i = 0; i < size; i++){
         Console.Write($"{data[i]:F12}");
         if (i < size - 1) Console.Write(", ");
      }
      Console.WriteLine("]");
   }
}

public class matrix{
   public readonly int size1, size2;
   private double[] data;
   public matrix(int n, int m){
      size1 = n; size2 = m;
      data = new double[size1 * size2];
   }
   public double this[int i, int j]{
      get => data[i + j * size1];
      set => data[i + j * size1] = value;
   }
   
   public void Print(string name = "Matrix"){
      for (int i = 0; i < size1; i++){
         for (int j = 0; j < size2; j++) Console.Write($"{this[i, j]:F12} ");
         Console.WriteLine();
      }
   }
   public matrix Copy(){
      matrix copy = new matrix(size1, size2);
      for (int i = 0; i < size1; i++){
         for (int j = 0; j < size2; j++){
            copy[i, j] = this[i, j];
         }
      }
      return copy;
   }
   public static matrix Multiply(matrix A, matrix B){
      if (A.size2 != B.size1) throw new Exception("Matrix dimensions do not match for multiplication");
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
}

public class QR{
   public matrix Q, R;

   public QR(matrix A){
      int n = A.size1, m = A.size2;
      Q = A.Copy();
      R = new matrix(m, m);
      for (int i = 0; i < m; i++){
         for (int j = 0; j < i; j++){
            double dot = 0;
            for (int k = 0; k < n; k++) dot += Q[k, j] * Q[k, i];
            R[j, i] = dot;
            for (int k = 0; k < n; k++) Q[k, i] -= dot * Q[k, j];
         }
         double norm = 0;
         for (int k = 0; k < n; k++) norm += Q[k, i] * Q[k, i];
         R[i, i] = Math.Sqrt(norm);
         for (int k = 0; k < n; k++) Q[k, i] /= R[i, i];
      }
   }

   public vector solve(vector b, matrix Qt){
      int n = Q.size1; // 4
      int m = R.size1; // 3
      if (b.size != n) throw new Exception("Vector b size does not match matrix dimensions");
      
      vector y = new vector(m);
      for (int i = 0; i < m; i++){
         double sum = 0;
         for (int k = 0; k < n; k++) sum += Qt[i, k] * b[k];  
         y[i] = sum;
      }

      y.Print("y");

      vector x = new vector(m);
      for (int i = m - 1; i >= 0; i--){
         double sum = (double)y[i];
         for (int k = i + 1; k < m; k++){
            sum -= (double)R[i, k] * (double)x[k];
         }
         if (Math.Abs(R[i, i]) < 1e-14) throw new Exception($"Singular matrix at row {i}: R[i, i] too small.");
         x[i] = sum / (double)R[i, i];
      }

      return x;
   }

   public matrix Inverse(){
      int n = Q.size1;
      matrix I = new matrix(n, n);
      matrix B = new matrix(n, n);
      
      for (int i = 0; i < n; i++) I[i, i] = 1.0;

      matrix Qt = new matrix(n, n);
      for (int i = 0; i < n; i++)
         for (int j = 0; j < n; j++)
            Qt[i, j] = Q[j, i];

      for (int i = 0; i < n; i++){
         vector ei = new vector(n);
         for (int j = 0; j < n; j++) ei[j] = I[j, i];
         vector x = solve(ei, Qt);
         for (int j = 0; j < n; j++) B[j, i] = x[j];
      }
      
      return B;
   }

   public double det(){
      double determinant = 1.0;
      for (int i = 0; i < R.size1; i++){
         determinant *= R[i, i];
      }
      return determinant;
   }
}

class Program{
   static void Main(string[] args){
      int maxN = 200;
      using (System.IO.StreamWriter file = new System.IO.StreamWriter("out.times.data")){
         for (int N = 100; N <= maxN; N += 20){
            Stopwatch stopwatch = Stopwatch.StartNew();
            int n = 4, m = 4;
            matrix A = new matrix(n, m);
            vector b = new vector(n);
            var rnd = new System.Random(1); 
            for (int i = 0; i < n; i++){
               for (int j = 0; j < m; j++) A[i, j] = rnd.NextDouble();
            }
            for (int i = 0; i < n; i++) b[i] = rnd.NextDouble();
            
            Console.WriteLine("Matrix A");
            A.Print();
           
            QR qr = new QR(A);
            
            Console.WriteLine("Matrix Q");
            qr.Q.Print("Matrix Q");
            
            Console.WriteLine("Matrix R");
            qr.R.Print();
            
            Console.WriteLine("QR = A");
            matrix QR_A = matrix.Multiply(qr.Q, qr.R);
            QR_A.Print();
            
            Console.WriteLine("Q^T * Q = I ");
            matrix QT = new matrix(m, n);
            for (int i = 0; i < m; i++)
               for (int j = 0; j < n; j++)
                  QT[i, j] = qr.Q[j, i];
            matrix QTQ = matrix.Multiply(QT, qr.Q);
            QTQ.Print();
            
            vector x = qr.solve(b, QT);
            
            x.Print("x");
            
            vector Ax = new vector(A.size1);
            
            for (int i = 0; i < A.size1; i++){
               double sum = 0;
               for (int j = 0; j < A.size2; j++){
                  sum += A[i, j] * x[j];
               }
               Ax[i] = sum;
            }
            
            Ax.Print("Ax");
            
            b.Print("b");
            
            Console.WriteLine("Inverse of A using QR Factorization:");
            matrix A_inv = qr.Inverse();
            A_inv.Print();
            
            Console.WriteLine("Checking that A * A_inv = I:");
            matrix IdentityCheck = matrix.Multiply(A, A_inv);
            IdentityCheck.Print();

            Console.WriteLine("Determinant of A:");
            Console.WriteLine(qr.det());

            stopwatch.Stop();
            
            double elapsed = stopwatch.Elapsed.TotalSeconds;
            file.WriteLine($"{N} {elapsed}");
            Console.WriteLine($"N={N}, Time={elapsed:F6} sec");
         }
      }
   }
}
