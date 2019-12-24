#include"header.h"

double ExactSolution(double x) {
	return 3 * sin(x) - 2;
}
double SecondDerative(double x, double y, double y1) {
	return (-3) * y1 - y + 9 * cos(x) - 2;
};
double FirstDerative(double x, double y, double y1) {
	return y1;
}
double Runge(double h) {
	double y = -2.0;
	double derative0 = 3.0;
	double k1 = 0, k2 = 0, k3 = 0, k4 = 0, l1 = 0, l2 = 0, l3 = 0, l4 = 0;
	for (double x = a; x < b; x += h) {
		k1 = h * SecondDerative(x, y, derative0);
		l1 = h * FirstDerative(x, y,  derative0);
		k2 = h * SecondDerative(x + h / 2.0, y + l1 / 2.0, derative0 + k1/2.0);
		l2 = h * FirstDerative(x + h / 2.0, y + l1 / 2.0, derative0 + k1 / 2.0);

		k3 = h * SecondDerative(x + h / 2.0, y + l2/ 2.0, derative0 + k2 / 2.0);
		l3 = h * FirstDerative(x + h / 2.0, y + l2 / 2.0, derative0 + k2 / 2.0);

		k4 = h * SecondDerative(x + h, y + l3, derative0 + k3);
		l4 = h * FirstDerative(x + h, y + l3, derative0 + k3);
		y = y + 1 / 6.0 * (l1 + 2 * l2 + 2 * l3 + l4);
		derative0 = derative0 + 1 / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
	}
	return y;
}

void CurveValue(double* h, double eps, double *delta) {
	double exact = ExactSolution(b);
	*delta = fabs(exact - Runge(*h));
	while (*delta > eps) {
		*h /= 2;
		*delta = fabs(exact - Runge(*h));
	}
}

double DoubleConversionMethod(double eps) {
	double h = sqrt(sqrt(eps));
	double r = 4.0;
	double yh = Runge(h);
	h /= 2;
	double yh2 = Runge(h);
	while (fabs(yh - yh2) / (pow(2, r) - 1) >  eps) {
		yh = yh2;
		h /= 2;
		yh2 = Runge(h);
	}
	return h;
}

double Runge(double h, double end) {
	double y = -2.0;
	double derative0 = 3.0;
	double k1 = 0, k2 = 0, k3 = 0, k4 = 0, l1 = 0, l2 = 0, l3 = 0, l4 = 0;
	for (double x = a; x < end; x += h) {
		k1 = h * SecondDerative(x, y, derative0);
		l1 = h * FirstDerative(x, y, derative0);
		k2 = h * SecondDerative(x + h / 2.0, y + l1 / 2.0, derative0 + k1 / 2.0);
		l2 = h * FirstDerative(x + h / 2.0, y + l1 / 2.0, derative0 + k1 / 2.0);

		k3 = h * SecondDerative(x + h / 2.0, y + l2 / 2.0, derative0 + k2 / 2.0);
		l3 = h * FirstDerative(x + h / 2.0, y + l2 / 2.0, derative0 + k2 / 2.0);

		k4 = h * SecondDerative(x + h, y + l3, derative0 + k3);
		l4 = h * FirstDerative(x + h, y + l3, derative0 + k3);
		y = y + 1 / 6.0 * (l1 + 2 * l2 + 2 * l3 + l4);
		derative0 = derative0 + 1 / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
	}
	return y;
}
