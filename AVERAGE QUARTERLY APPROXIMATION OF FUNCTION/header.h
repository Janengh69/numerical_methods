#pragma once
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#define a 1
#define b 10

#define eps 1e-2
#define maxn 100
using namespace std;

double func(double x);

void PrintMatrix(double** matrix, int n);
double Req(double x, int n);

double Pm(double x, int m, double* matrix_roots);

double PowReq(double x, int i, int j);

double SimpsoMethod(double h, int n, int i, int j);

double Integral(int i, int j);

double** CreateMatrix(int n);
double getDelta(int m, double* matrix_roots);
double* GaussElimination(double** initial_matrix, int n);
