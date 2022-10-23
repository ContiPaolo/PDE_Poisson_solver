# include <iostream>
using namespace std;
#include <utility>
void validation1_create(double * & anal, double * & b, int N, int M);
void validation2_create(double * & anal, double * & b, int N, int M);
void validation3_create(double * & anal, double * & b, int N, int M);
void validation4_create(double * & anal, double * & b, int N, int M);
void validation5_create(double * & anal, double * & b, int N, int M);
void validation6_create(double * & anal, double * & b, int N, int M);
double errl2(double * & x, double *& anal, int N, int M);
double errlinf(double * & x, double *& anal, int N, int M);
pair<int,int> bij(int k, int N, int M);
