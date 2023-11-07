extern "C" {
#include "MatrixIO.h"
#include "MatrixUtils.h"
}
#include <cstdlib>
#include <random>
#include <iostream>
#include <algorithm>
#include <set>
using namespace std;

void generate_basic_square_sparce_matrix(int size, const char* filename_out) {
	const int k = min(25, size);
	matrix_COO mtx;
	mtx.N = size;
	mtx.M = mtx.N;
	mtx.nz = 25 * mtx.N;
	mtx.I = (int*)malloc(mtx.nz * sizeof(int));
	mtx.J = (int*)malloc(mtx.nz * sizeof(int));
	mtx.val = (double*)malloc(mtx.nz * sizeof(double));

	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<int> rand_col(0, mtx.M - 1);
	uniform_real_distribution<double> rand_val(-9.99, 9.99);

	for (int i = 0; i < mtx.N; i++) {
		int j = i * k;
		set<int> was;
		for (int y = 0; y < k; y++) {
			int col;
			do {
				col = rand_col(gen);
			} while (was.count(col));
			mtx.I[j + y] = i;
			mtx.J[j + y] = col;
			mtx.val[j + y] = rand_val(gen);
		}
	}
	save_matrix_BIN(filename_out, &mtx);
	delete_COO(&mtx);
}

int main() {
	generate_basic_square_sparce_matrix(20'000, "D:\\source\\coursework\\sample_matrices\\multA.bin");
	generate_basic_square_sparce_matrix(20'000, "D:\\source\\coursework\\sample_matrices\\multB.bin");
	return 0;
}
