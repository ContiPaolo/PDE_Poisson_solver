#include "Validation.hpp"
#include <iostream>
using namespace std;
#include <utility>
#include <cmath>
#include <algorithm>

#define PI 3.14159265

//First test function: phi(x,y)=sin(PI*x)*sin(PI*y)
void validation1_create(double * & anal, double * & b,  int N, int M){
  double h1 = 1./N;
  double h2 = 1./M;
  for(int k = 0; k<(N-1)*(M-1); ++k){
    pair<int,int> v = bij(k, N, M);
    double n = get<0>(v); double x_n = h1*n; //n = 1 : N-1
    double m = get<1>(v); double y_m = h2*m; //m = 1 : M-1
    anal[k] = sin(PI*x_n)*sin(PI*y_m);
    b[k] = 2.*PI*PI*sin(PI*x_n)*sin(PI*y_m);
  }
}

//Second test function: phi(x,y)=xy(x-1)(y-1)
void validation2_create(double * & anal, double * & b,  int N, int M){
  double h1 = 1./N;
  double h2 = 1./M;
  for(int k = 0; k<(N-1)*(M-1); ++k){
    pair<int,int> v = bij(k, N, M);
    double n = get<0>(v); double x_n = h1*n; //n = 1 : N-1
    double m = get<1>(v); double y_m = h2*m; //m = 1 : M-1
    anal[k] = x_n*y_m*(1.-y_m)*(1.-x_n);
    b[k] = 2.*(y_m*(1.-y_m) + x_n*(1.-x_n));
  }
}

//Third test function: phi(x,y) = exp(-1000*pow(x-0.5,2)-1000*pow(y-0.5,2))
void validation3_create(double * & anal, double * & b,  int N, int M){
  double h1 = 1./N;
  double h2 = 1./M;
  for(int k = 0; k<(N-1)*(M-1); ++k){
    pair<int,int> v = bij(k, N, M);
    double n = get<0>(v); double x_n = h1*n; //n = 1 : N-1
    double m = get<1>(v); double y_m = h2*m; //m = 1 : M-1
    anal[k] = exp(-1000*pow(x_n-0.5,2)-1000*pow(y_m-0.5,2));
    b[k] = -exp(-1000*pow(x_n-0.5,2)-1000*pow(y_m-0.5,2))*(4*1e6*(x_n*x_n-x_n+y_m*y_m-y_m) + 1.996*1e6);
  }
}

//Fourth test function: phi(x,y)= x*y^0.01*(1-y)*(1-x)
void validation4_create(double * & anal, double * & b,  int N, int M){
  double h1 = 1./N;
  double h2 = 1./M;
  for(int k = 0; k<(N-1)*(M-1); ++k){
    pair<int,int> v = bij(k, N, M);
    double n = get<0>(v); double x_n = h1*n; //n = 1 : N-1
    double m = get<1>(v); double y_m = h2*m; //m = 1 : M-1
    anal[k] = x_n*pow(y_m,0.6)*(1.-y_m)*(1.-x_n); //x*y^0.01*(1-y)*(1-x)
    b[k] = +(0.24*(x_n-1.)*x_n*(y_m-1.))/pow(y_m,1.4) - (1.2*(x_n-1.)*x_n)/pow(y_m,0.4) - 2*(y_m-1.)*pow(y_m,0.6);
  }
}

//Fifth test function: phi(x,y)=x*y^0.2*(1-y)*(1-x)
void validation5_create(double * & anal, double * & b,  int N, int M){
  double h1 = 1./N;
  double h2 = 1./M;
  for(int k = 0; k<(N-1)*(M-1); ++k){
    pair<int,int> v = bij(k, N, M);
    double n = get<0>(v); double x_n = h1*n; //n = 1 : N-1
    double m = get<1>(v); double y_m = h2*m; //m = 1 : M-1
    anal[k] = x_n*pow(y_m,0.2)*(1.-y_m)*(1.-x_n);
    b[k] = +(0.16*(x_n-1.)*x_n*(y_m-1.))/pow(y_m,1.8) - (0.4*(x_n-1.)*x_n)/pow(y_m,0.8) - 2*(y_m-1.)*pow(y_m,0.2);
  }
}

//Sixth test function: phi(x,y)=10*x*y^0.01*(1-y)*(1-x)
void validation6_create(double * & anal, double * & b,  int N, int M){
  double h1 = 1./N;
  double h2 = 1./M;
  for(int k = 0; k<(N-1)*(M-1); ++k){
    pair<int,int> v = bij(k, N, M);
    double n = get<0>(v); double x_n = h1*n; //n = 1 : N-1
    double m = get<1>(v); double y_m = h2*m; //m = 1 : M-1
    anal[k] = 10*x_n*pow(y_m,0.01)*(1.-y_m)*(1.-x_n); //10*x*y^0.01*(1-y)*(1-x)
    b[k] = +(0.099*(x_n-1.)*x_n*(y_m-1.))/pow(y_m,1.99) - (0.2*(x_n-1.)*x_n)/pow(y_m,0.99) - 20*(y_m-1.)*pow(y_m,0.01);
  }
}

//the "l" function as defined in our latex
pair<int,int> bij(int k, int N, int M){
  int x = k+1 - int(k/(N-1))*(N-1);
  int y = M-1 - int(k/(N-1));
  return make_pair(x,y);
}

//l2 error:
double errl2(double * &x, double *& anal, int N, int M){
  double h1 = 1./N; double h2 = 1./M;
  double err = 0.;
  for(int k = 0; k<(N-1)*(M-1); ++k){
    err = err+ pow(anal[k]-x[k],2.);
  }
    return sqrt(h1*h2*err);
}

//l-infinity error:
double errlinf(double * &x, double *& anal, int N, int M){
  double h1 = 1./N; double h2 = 1./M;
  double m = 0.;
  for(int k = 0; k<(N-1)*(M-1); ++k){
    m = max(abs(anal[k]-x[k]), m);
  }
  return m;
}
