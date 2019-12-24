#include<iostream>
#include<cmath>
#include <iomanip>

#define a 0 
#define M_PI 3.14159265358979323846
#define b M_PI/2

using namespace std;

double ExactSolution(double x);
double SecondDerative(double x, double y, double first_derative);
double FirstDerative(double x, double y, double y1);
double Runge(double h); 
double DoubleConversionMethod(double eps);
double Runge(double h, double end);
void CurveValue(double* h, double eps, double*delta);