#include <iostream>
extern "C" {
	#include "mmio.h"
}
using namespace std;

int main() {
	const char* filename = "D:\\source\\coursework\\ASIC_320k.mtx";
	int M, N, nz;
	int* I;
	int* J;
	double* val;
	if (mm_read_unsymmetric_sparse(filename, &M, &N, &nz, &val, &I, &J) != 0) {
		cout << "error in reading matrix" << endl;
		return -1;
	}
	cout << "matrix read succsessfully" << endl;
	cout << "N: " << N << ", M: " << M << ", nz:" << nz << endl;
	double* V = new double[M];

	double* Result = new double[N];

	delete[] I;
	delete[] J;
	delete[] val;
	delete[] V;
	delete[] Result;
	return 0;
}
