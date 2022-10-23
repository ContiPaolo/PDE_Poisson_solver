#include <iostream>
#include "Stencils.hpp"
#include "Validation.hpp"
#include "PlotFunction.hpp"
//#include "/opt/local/include/umfpack.h" //for Yvonne
#include <suitesparse/umfpack.h> //for Paolo
#include <stdio.h>
#include <cassert>
#include <fstream>
#include <string>
#include <chrono>
using namespace std;
using namespace chrono;

int main(){
  void *Symbolic, *Numeric;
  static const int calls = 8; //number of times we will test our code with increasing number of discretisations in [0,1].
  int N[calls] = {5,  7,  11,  18, 27, 42, 65, 100};// Chose the numbers spaced evenly on a log scale
  int M[calls] = {5,  7,  11,  18, 27, 42, 65, 100}; // Take the number of discretisations in x and y equal for the analysis
  //of convergence. Otherwise, the code works for N and M possibily different.
  double e[calls];
  int duration[calls];

  for(int k=0; k<calls; k++){
    auto start = high_resolution_clock::now(); //start to measure time
    Stencils A(N[k],M[k]);
    int dim = A.m_dim; // the matrix we are inverting, A,  belongs in dim by dim.
    int nz = A.num_elem; // total number of non-zero elements
    double* b = new double[dim]; double * x = new double[dim]; double * anal = new double[dim];

    //populate the analytical and RHS vector:
    validation1_create(anal, b, N[k], M[k]);
    printf("%i: The matrix is %d x %d , total number of non zero elements is %d \n", k, dim, dim, nz);

    //Allocating memory for the CCS vectors:
    double * Ax = new double[nz];
    int * Ap = new int[dim+1]; Ap[0]=0;
    int * Ai = new int[nz];
    CCS(A, Ax, Ap, Ai); //populates the CCS vectors

    umfpack_di_symbolic(dim, dim, Ap, Ai, Ax, &Symbolic, NULL, NULL);
    umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL);
    umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, NULL, NULL); // computed solution placed in x

    e[k] = errl2(x,anal, N[k], M[k]);
    //To visualize the error:
    //printf("%i: l2 error is %f \n", k, e[k]);
    //e[k] = errlinf(x,anal, N[k], M[k]);

    umfpack_di_free_symbolic(&Symbolic);
    umfpack_di_free_numeric(&Numeric);

    delete[] Ax; delete[] Ap; delete[] Ai; delete[] b; //free allocated memory

    auto end = high_resolution_clock::now(); //end of mesuring time
    duration[k] = duration_cast<microseconds> (end-start).count();
    printf("%i: Duration: %i ns \n", k, duration[k]);

    if(k==calls-1){ //I create the file to plot the functions only for the last k
    char namefile[20];
    sprintf(namefile, "Fonction1_%d.txt", k);
    PlotFunction(namefile, anal, x, N[k], M[k]);
    }

    delete[] x;
  }

  //I create the files to plot the error and duration:
  PlotError("Error1.txt", e, N, M, calls);
  PlotDuration("Duration1.txt", duration, N, M, calls);

  return 0;
}
