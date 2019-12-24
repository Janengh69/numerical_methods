
#include"header.h"

int main() {
	double h = (b - a) / 10;
	cout << "Differential equation: y'' + 3*y' + y - 9 * cos(x) + 2 = 0" << endl << endl;
	double eps = 1e-2;
	double delta = 0;
	CurveValue(&h, eps, &delta);
	cout << "Delta " << delta << " for Runge-Kutta " << Runge(h) << " with h: " << h;
	cout << "\n\nH with the same error = " << DoubleConversionMethod(delta) << endl;
	h = (b - a) / 10;
	cout << "\nTable of values of the integral curve: " << endl;
	cout << "\   X\t   Y\n";
	for (double x = a; x <= b; x += h) {
		cout << setprecision(4) << x << "\t" << setprecision(6) << Runge(h, x) << endl;
	}
	return 0;
}