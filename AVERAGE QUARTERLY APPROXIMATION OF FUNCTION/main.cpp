#include "header.h"

int main() {

	double** matrix;
	double delta;
	double* matrix_roots;

	ofstream out_file("graph.csv");

	matrix = CreateMatrix(maxn);
	for (int i = 1; i < maxn; i++) {
		matrix_roots = GaussElimination(matrix, i);
		delta = getDelta(i + 1, matrix_roots);
		cout << "Delta = " << delta << " m = " << i + 1  << endl;
		if (delta < eps) {
			double h = (b - a) / 100.0;
			double x = a;
			cout << "\nLeast Squares Deviation = " << delta << " for m = " << i << endl;
			for (int j = 0; j <= maxn; j++) {
				x = a + j * h;
				cout << "X = "<< setprecision(2) << x << "  Pm  = " << setprecision(3) << Pm(x, i + 1, matrix_roots) << endl;
				out_file << setprecision(3) << x << "," << setprecision(5) << Pm(x, i + 1, matrix_roots) << endl;
			}
			break;
		}
	}
	system("pause");
	return 0;
}