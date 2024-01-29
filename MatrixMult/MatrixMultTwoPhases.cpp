#include "MatrixMult.h"
extern "C" {
#include "MatrixUtils.h"
}
#include <stdlib.h>
#include <vector>
#include <queue>
#include <utility>
#include <algorithm>

int matrix_mult_TF_first_phase(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr, int delta) {
	int N = A_csr->N;
	std::vector<std::vector<int>> col(N);

	std::queue<std::pair<int, int>> tasks;
	for (int i = 0; i < N; i += delta) {
		tasks.push({ i, std::min(i + delta, N) });
	}

#pragma omp parallel
	{
		int A = 0, B = 1;
		std::vector<int> X(N, -1);
		while (A < B) {
#pragma omp critical
			{
				if (!tasks.empty()) {
					auto t = tasks.front(); tasks.pop();
					A = t.first;
					B = t.second;
				}
				else {
					A = B = 0;
				}
			}
			for (int i = A; i < B; i++) {
				for (int j = 0; j < N; j++) {
					int temp = -1;
					for (int k = A_csr->row_id[i]; k < A_csr->row_id[i + 1]; k++) {
						X[A_csr->col[k]] = k;
					}
					for (int l = B_csr->row_id[j]; temp == -1 && l < B_csr->row_id[j + 1]; l++) {
						temp = X[B_csr->col[l]];
					}
					for (int k = A_csr->row_id[i]; k < A_csr->row_id[i + 1]; k++) {
						X[A_csr->col[k]] = -1;
					}
					if (temp != -1) {
						col[i].push_back(j);
					}
				}
			}
		}
	}
	Res_csr->N = N;
	Res_csr->M = N;
	Res_csr->row_id = (int*)malloc((N + 1) * sizeof(int));
	int nz = 0;
	for (int i = 0; i < N; i++) {
		Res_csr->row_id[i] = nz;
		nz += col[i].size();
	}
	Res_csr->row_id[N] = nz;
	Res_csr->col = (int*)malloc(nz * sizeof(int));
	Res_csr->value = (double*)malloc(nz * sizeof(double));
	for (int i = 0; i < N; i++) {
		if (col[i].size()) {
			std::memcpy(Res_csr->col + Res_csr->row_id[i], col[i].data(), sizeof(int) * col[i].size());
		}
	}
	return 0;
}

int matrix_mult_TF_second_phase(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr, int delta) {
	int N = A_csr->N;

	std::queue<std::pair<int, int>> tasks;
	for (int i = 0; i < N; i += delta) {
		tasks.push({ i, std::min(i + delta, N) });
	}

#pragma omp parallel
	{
		int A = 0, B = 1;
		std::vector<int> X(N, -1);
		while (A < B) {
#pragma omp critical
			{
				if (!tasks.empty()) {
					auto t = tasks.front(); tasks.pop();
					A = t.first;
					B = t.second;
				}
				else {
					A = B = 0;
				}
			}
			for (int i = A; i < B; i++) {
				for (int u = Res_csr->row_id[i]; u < Res_csr->row_id[i + 1]; u++) {
					int j = Res_csr->col[u];
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
					Res_csr->value[u] = t;
				}
			}
		}
	}
	return 0;
}
