#include<fstream>
#include<iostream>
#include<cassert>
#include<string>

void PlotFunction(char * namefile, double * & anal, double * & x,int N, int M){
  ofstream file;
  file.open(namefile);
  assert(file.is_open());
  double h1 = 1./N;
  double h2 = 1./M;

  //plot of the inner domain
  for(int k = 0; k<(N-1)*(M-1); ++k){
    pair<int,int> v = bij(k, N, M);
    double n = get<0>(v); double x_n = h1*n; //n = 1 : N-1
    double m = get<1>(v); double y_m = h2*m; //m = 1 : M-1
    file << x_n << '\t' << y_m << '\t' << anal[k] << '\t' << x[k];
    file << std::endl << std::endl;
  }

  //plot of the boundary
  for(int k=0; k*h1<=1; k++){
    file << k*h1 << '\t' << 0 << '\t' << 0 << '\t' << 0;
    file << std::endl << std::endl;
    file << k*h1 << '\t' << 1 << '\t' << 0 << '\t' << 0;
    file << std::endl << std::endl;
  }
  for(int k=h2; k*h2<1; k++){
    file << 0 << '\t' << k*h1 << '\t' << 0 << '\t' << 0;
    file << std::endl << std::endl;
    file << 1 << '\t' << k*h2 << '\t' << 0 << '\t' << 0;
    file << std::endl << std::endl;
  }
  file.close();
}

void PlotError(string namefile, double * err, int * N, int * M, const int calls){
  ofstream file;
  file.open(namefile);
  assert(file.is_open());
  for(int k=0; k<calls; k++){
    assert( N[k] == M[k]);
    double h = 1./N[k];
    file << h  << '\t' << h*h  << '\t' << err[k];
    file <<std::endl;
  }
  file.close();
}

void PlotDuration(string namefile, int * duration, int * N, int * M, const int calls){
  ofstream file;
  file.open(namefile);
  assert(file.is_open());
  for(int k=0; k<calls; k++){
    assert(N[k] == M[k]);
    file << N[k]  << '\t' << N[k]*N[k] << '\t' << pow(N[k],3)  << '\t' << duration[k];
    file <<std::endl;
  }
  file.close();
}
/* GNUPLOT:
#### FUNCTIONS: #####
set xrange[0:1]
set xlabel "x"
set ylabel "y"
set yrange[0:1]
set ticslevel 0
splot "Fonction1_7.txt" using 1:2:3 with lines title "analytical solution",\
"Fonction1_7.txt" using 1:2:4 with lines title "numerical solution"

#### ERRORS: #####
set logscale xy
set xlabel "h"
set key left top
plot "Error1.txt" u 1:1 w lp title "h",\
"Error1.txt" u 1:2 w lp title "h^2",\
"Error1.txt" u 1:3 w lp title "l2 error"

#### DURATION: #####
set logscale xy
set xlabel "N"
set key left top
set ylabel "nanoseconds"
plot "Duration2.txt" u 1:1 w lp title "N",\
"Duration1.txt" u 1:2 w lp title "N^2",\
"Duration1.txt" u 1:3 w lp title "N^3",\
"Duration1.txt" u 1:4 w lp title "Duration"
*/
