using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Lab6
{
    class Program
    {
        static void Main(string[] args)
        {
            IntCubicSpline s = new IntCubicSpline(2, 10,50);
            s.Out();
        }
    }
}
