using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace Lab6
{
    public class NoRootsException : System.ApplicationException
    {
        public NoRootsException() { }
        public NoRootsException(string message) : base(message) { }
    }

    public class OutOfIntervalException : System.ApplicationException
    {
        public OutOfIntervalException() { }
        public OutOfIntervalException(string message) : base(message) { }
    }

    class IntCubicSpline
    {
        static int N;
        double x1, x2,h;
        double[] B;
        double[,] A;
        double[] X;
        delegate double Funct(double x);
        
        public IntCubicSpline(double X1, double X2,int P)
        {
            x1 = X1;
            x2 = X2;
            N = P;
            B = new double[N];
            A = new double [N,N];
            X = new double[N];//Ci
            h = (x2 - x1) / (double)N;
        }

        double F(double x)//function for interpolation
        {
            return Math.Log10(x) * Math.Log(10 * x) * Math.Sin(2.5 * x);
        }

        void PrintMatrix(int n,double[,] Ac, double[] Bc)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Console.Write("{0,5}", Ac[i, j]);
                }
                Console.Write(" |{0,5}", Bc[i]);
                Console.WriteLine();
            }
        }

        void PrintRoots(int n, double[] Bc)
        {
            for (int i = 0; i < n; i++)
            {
               Console.Write("x{0} = {1,10}", i + 1, Bc[i]);
               Console.WriteLine();
            }
        }

        void GElimination(int n, double[,] _Ac, double[] _Bc)//Gaussian elimination algorythm for matrix roots searching
        {
            double[,] Ac = new double[n,n];
            double[] Bc = new double[n];
            Array.Copy(_Ac, Ac, _Ac.LongLength);
            _Bc.CopyTo(Bc,0);
            double[] Buf = new double[n+1];
            for (int k = 0; k < n; k++)
            {
                if (Ac[k, k] == 0) { throw new NoRootsException("There is no roots in this system of equations !"); }
                else
                {
                    double elem = Ac[k, k];
                    for (int j = 0; j < n; j++)
                    {
                        Ac[k, j] /= elem;
                        Buf[j] = Ac[k, j];
                    }
                    Bc[k] /= elem;
                    Buf[n] = Bc[k];
                    for (int i = k + 1; i < n;i++ )
                    {
                        double coef = Ac[i, k];
                        for (int j = 0; j < n; j++)
                        {
                            Ac[i, j] -= Buf[j] * coef;
                        }
                        Bc[i] -= Buf[n] * coef;
                    }
                }
            }
            BackSub(n, Ac, Bc);
        }

        void BackSub(int n, double[,] Ac, double[] Bc)//Backward substitution for previous method
        {
            double summ;
            for (int k = n - 1; k >= 0;k--)
            {
                summ = 0;
                for (int j = k + 1; j < n; j++)
                {
                    summ += Ac[k, j] * Bc[j];
                }
                Bc[k] = Bc[k] - summ;
            }
            X = Bc;
        }


        double Ai(double a,double b,int i,Funct f)
        {
            return f(a+i*h);
        }

        void FillM(int n, double a, double b,Funct f)
        {
            A[0, 0] = 1;
            A[n-1, n-1] = 1;
            for (int i = 1; i < n-1; i++)
            {
                A[i,i]=4*h;
                A[i, i - 1] = h;
                if (i < n - 1)
                    A[i, i + 1] = h;
                double fi0 = f(a+h*(i-1));
                double fi = f(a+h*i);
                double fi1 = f(a+h*(i+1));
                double s1 = (fi1+fi)/h;
                double s2 = (fi-fi0)/h;
                B[i]=6*(s1-s2);
            }
        }

        double Di(int i)
        {
            if (i == 0)
            {
                return X[i] / h;
            }
            return (X[i] - X[i - 1]) / h;
        }

        double Bi(int i,double a,Funct f)
        {
            double fi0 = f(a + h * (i - 1));
            double fi = f(a + h * i);
            double di = Di(i);
            return (h * X[i]) / 2 - (di * h * h) / 6 + (fi - fi0) / h;           
        }

        double Si(double x,double a, double b,Funct f)
        {
            if ((x < a) || (x > b)) { throw new OutOfIntervalException("X is out of range of [x1,x2]!"); }
            if (x == a) { return Ai(a,b, 0, f); }
            int i = 0;
            double xi;
            do
            {
                i++;
                xi = a+i*h;
            }while (x>xi);
            double ai = Ai(a,b,i,f);
            double bi = Bi(i,a,f);
            double ci = X[i];
            double di = Di(i);
            double delta = x-xi;
            return ai + bi * delta + (ci * delta * delta) / 2 + (di*delta*delta*delta)/2;
        }

        void PrintSpline(int n,double a, double b,Funct f)
        {
            double x = 0.2,y;
            StreamWriter sw = new StreamWriter("Spline.txt", false,Encoding.Default);
            for (int i = 0; i < n; i++)
            {
                x = a + i * h;
                y = Si(x, a, b, f);
                sw.WriteLine("{0,6} {1}", x, y);
                Console.WriteLine("{0};{1}",x,y);
            }
            sw.Close();
        }

        public void Out()
        {          
            FillM(N, x1, x2, F);
            GElimination(N, A, B);
            PrintSpline(N, x1, x2, F);
        }


    }
}
