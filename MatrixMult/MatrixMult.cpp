#include "MatrixMult.h"
#include <stdlib.h>
#include <vector>

int matrix_mult_naive_1(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr) {
	// A[N][N] * B[N][N] = C[N][N]
	int N = A_csr->N;
	int nz = 0;
	std::vector<int> row;
	std::vector<int> col;
	std::vector<double> val;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double t = 0.0;
			for (int k = A_csr->row_id[i]; k < A_csr->row_id[i + 1]; k++) {
				for (int l = B_csr->row_id[j]; l < B_csr->row_id[j + 1]; l++) {
					if (A_csr->col[k] == B_csr->col[l]) {
						t += A_csr->value[k] * B_csr->value[l];
					}
				}
			}
			if (t != 0) {
				row.push_back(i);
				col.push_back(j);
				val.push_back(t);
				nz++;
			}
		}
	}
	matrix_COO res_coo;
	res_coo.N = res_coo.M = N;
	res_coo.nz = nz;
	res_coo.I = row.data();
	res_coo.J = col.data();
	res_coo.val = val.data();
	convert_COO_to_CSR(&res_coo, Res_csr);
	return 0;
}

int matrix_mult_naive_2(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr) {
	// A[N][N] * B[N][N] = C[N][N]
	int N = A_csr->N;
	int nz = 0;
	std::vector<int> row;
	std::vector<int> col;
	std::vector<double> val;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double t = 0.0;
			for (int k = A_csr->row_id[i], l = B_csr->row_id[j]; k < A_csr->row_id[i + 1] && l < B_csr->row_id[j + 1];) {
				if (A_csr->col[k] < B_csr->col[l]) {
					k++;
				}
				else if (A_csr->col[k] > B_csr->col[l]) {
					l++;
				}
				else { // A_csr->col[k] == B_csr->col[l]
					t += A_csr->value[k] * B_csr->value[l];
					k++;
					l++;
				}
			}
			if (t != 0) {
				row.push_back(i);
				col.push_back(j);
				val.push_back(t);
				nz++;
			}
		}
	}
	matrix_COO res_coo;
	res_coo.N = res_coo.M = N;
	res_coo.nz = nz;
	res_coo.I = row.data();
	res_coo.J = col.data();
	res_coo.val = val.data();
	convert_COO_to_CSR(&res_coo, Res_csr);
	return 0;
}

int matrix_mult_naive_3(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr) {
	// A[N][N] * B[N][N] = C[N][N]
	int N = A_csr->N;
	int nz = 0;
	std::vector<int> row;
	std::vector<int> col;
	std::vector<double> val;
	std::vector<int> X(N, -1);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double t = 0.0;
			for (int k = A_csr->row_id[i]; k < A_csr->row_id[i + 1]; k++) {
				X[A_csr->col[k]] = k;
			}
			for (int l = B_csr->row_id[j]; l < B_csr->row_id[j + 1]; l++) {
				if (X[B_csr->col[l]] != -1) {
					t += A_csr->value[X[B_csr->col[l]]] * B_csr->value[l];
				}
			}
			for (int k = A_csr->row_id[i]; k < A_csr->row_id[i + 1]; k++) {
				X[A_csr->col[k]] = -1;
			}
			if (t != 0) {
				row.push_back(i);
				col.push_back(j);
				val.push_back(t);
				nz++;
			}
		}
	}
	matrix_COO res_coo;
	res_coo.N = res_coo.M = N;
	res_coo.nz = nz;
	res_coo.I = row.data();
	res_coo.J = col.data();
	res_coo.val = val.data();
	convert_COO_to_CSR(&res_coo, Res_csr);
	return 0;
}

int matrix_mult_naive_1_naiveomp(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr) {
	// A[N][N] * B[N][N] = C[N][N]
	int N = A_csr->N;
	int nz = 0;
	std::vector<int> row;
	std::vector<int> col;
	std::vector<double> val;
#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double t = 0.0;
			for (int k = A_csr->row_id[i]; k < A_csr->row_id[i + 1]; k++) {
				for (int l = B_csr->row_id[j]; l < B_csr->row_id[j + 1]; l++) {
					if (A_csr->col[k] == B_csr->col[l]) {
						t += A_csr->value[k] * B_csr->value[l];
					}
				}
			}
			if (t != 0) {
#pragma omp critical
				{
					row.push_back(i);
					col.push_back(j);
					val.push_back(t);
					nz++;
				}
			}
		}
	}
	matrix_COO res_coo;
	res_coo.N = res_coo.M = N;
	res_coo.nz = nz;
	res_coo.I = row.data();
	res_coo.J = col.data();
	res_coo.val = val.data();
	convert_COO_to_CSR(&res_coo, Res_csr);
	return 0;
}

int matrix_mult_naive_2_naiveomp(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr) {
	// A[N][N] * B[N][N] = C[N][N]
	int N = A_csr->N;
	int nz = 0;
	std::vector<int> row;
	std::vector<int> col;
	std::vector<double> val;
#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double t = 0.0;
			for (int k = A_csr->row_id[i], l = B_csr->row_id[j]; k < A_csr->row_id[i + 1] && l < B_csr->row_id[j + 1];) {
				if (A_csr->col[k] < B_csr->col[l]) {
					k++;
				}
				else if (A_csr->col[k] > B_csr->col[l]) {
					l++;
				}
				else { // A_csr->col[k] == B_csr->col[l]
					t += A_csr->value[k] * B_csr->value[l];
					k++;
					l++;
				}
			}
			if (t != 0) {
#pragma omp critical
				{
					row.push_back(i);
					col.push_back(j);
					val.push_back(t);
					nz++;
				}
			}
		}
	}
	matrix_COO res_coo;
	res_coo.N = res_coo.M = N;
	res_coo.nz = nz;
	res_coo.I = row.data();
	res_coo.J = col.data();
	res_coo.val = val.data();
	convert_COO_to_CSR(&res_coo, Res_csr);
	return 0;
}

int matrix_mult_naive_3_naiveomp(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr) {
	// A[N][N] * B[N][N] = C[N][N]
	int N = A_csr->N;
	int nz = 0;
	std::vector<int> row;
	std::vector<int> col;
	std::vector<double> val;
	std::vector<int> X(N, -1);
#pragma omp parallel
	{
		std::vector<int> X(N, -1);
#pragma omp for
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				double t = 0.0;
				for (int k = A_csr->row_id[i]; k < A_csr->row_id[i + 1]; k++) {
					X[A_csr->col[k]] = k;
				}
				for (int l = B_csr->row_id[j]; l < B_csr->row_id[j + 1]; l++) {
					if (X[B_csr->col[l]] != -1) {
						t += A_csr->value[X[B_csr->col[l]]] * B_csr->value[l];
					}
				}
				for (int k = A_csr->row_id[i]; k < A_csr->row_id[i + 1]; k++) {
					X[A_csr->col[k]] = -1;
				}
				if (t != 0) {
#pragma omp critical
					{
						row.push_back(i);
						col.push_back(j);
						val.push_back(t);
						nz++;
					}
				}
			}
		}
	}
	matrix_COO res_coo;
	res_coo.N = res_coo.M = N;
	res_coo.nz = nz;
	res_coo.I = row.data();
	res_coo.J = col.data();
	res_coo.val = val.data();
	convert_COO_to_CSR(&res_coo, Res_csr);
	return 0;
}
