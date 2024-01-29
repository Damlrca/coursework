extern "C" {
#include "MatrixIO.h"
#include "MatrixUtils.h"
}
#include <cstdlib>
#include <random>
#include <iostream>
#include <algorithm>
#include <set>

#include "MatrixGenerator.hpp"

using namespace std;

random_device rd;
mt19937 gen(rd());

set<int> random_sample(int n, int m) {
	set<int> res;
	for (int r = n - m; r < n; r++) {
		int i = uniform_int_distribution<int>(0, r)(gen);
		if (res.count(i))
			res.insert(r);
		else
			res.insert(i);
	}
	return res;
}

void generate_uniform_square_sparse_matrix(int n, int k, matrix_CSR& res) {
	uniform_real_distribution<double> rand_val(-9.99, 9.99);
	matrix_COO mtx;
	mtx.N = n;
	mtx.M = mtx.N;
	mtx.nz = k * mtx.N;
	mtx.I = (int*)malloc(mtx.nz * sizeof(int));
	mtx.J = (int*)malloc(mtx.nz * sizeof(int));
	mtx.val = (double*)malloc(mtx.nz * sizeof(double));

	int j = 0;
	for (int i = 0; i < mtx.N; i++) {
		set<int> cols = random_sample(mtx.N, k);
		for (int col : cols) {
			mtx.I[j] = i;
			mtx.J[j] = col;
			mtx.val[j] = rand_val(gen);
			j++;
		}
	}

	convert_COO_to_CSR(&mtx, &res);
	delete_COO(&mtx);
}

void generate_nonuniform_square_sparse_matrix(int n, int k1, int k2, matrix_CSR& res) {
	uniform_real_distribution<double> rand_val(-9.99, 9.99);
	int nz = 0;
	for (int i = 0; i < n; i++) {
		int k = i * (k2 - k1) / (n - 1) + k1;
		nz += k;
	}

	matrix_COO mtx;
	mtx.N = n;
	mtx.M = mtx.N;
	mtx.nz = nz;
	mtx.I = (int*)malloc(mtx.nz * sizeof(int));
	mtx.J = (int*)malloc(mtx.nz * sizeof(int));
	mtx.val = (double*)malloc(mtx.nz * sizeof(double));

	int j = 0;
	for (int i = 0; i < mtx.N; i++) {
		int k = i * (k2 - k1) / (n - 1) + k1;
		set<int> cols = random_sample(mtx.N, k);
		for (int col : cols) {
			mtx.I[j] = i;
			mtx.J[j] = col;
			mtx.val[j] = rand_val(gen);
			j++;
		}
	}

	convert_COO_to_CSR(&mtx, &res);
	delete_COO(&mtx);
}

// -------------------------------------------------------------------------------------
// old function:
void generate_basic_square_sparce_matrix_BIN_FILE(int size, int k, const char* filename_out) {
	k = min(k, size);
	matrix_COO mtx;
	mtx.N = size;
	mtx.M = mtx.N;
	mtx.nz = k * mtx.N;
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
			was.insert(col);
			mtx.I[j + y] = i;
			mtx.J[j + y] = col;
			mtx.val[j + y] = rand_val(gen);
		}
	}
	save_matrix_BIN(filename_out, &mtx);
	cout << filename_out << " " << "generated and saved" << endl;
	delete_COO(&mtx);
}

/*
int main() {
	generate_basic_square_sparce_matrix_BIN_FILE(10, 2, "D:\\source\\coursework\\sample_matrices\\multA_small.bin");
	generate_basic_square_sparce_matrix_BIN_FILE(10, 2, "D:\\source\\coursework\\sample_matrices\\multB_small.bin");
	generate_basic_square_sparce_matrix_BIN_FILE(20'000, 25, "D:\\source\\coursework\\sample_matrices\\multA.bin");
	generate_basic_square_sparce_matrix_BIN_FILE(20'000, 25, "D:\\source\\coursework\\sample_matrices\\multB.bin");
	return 0;
}
*/
