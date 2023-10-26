#include "MatrixUtils.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void delete_matrix(int* N, int* M, int* nz, double** val, int** I, int** J) {
	*N = 0;
	*M = 0;
	*nz = 0;
	free(*val); *val = NULL;
	free(*I); *I = NULL;
	free(*J); *J = NULL;
}

void delete_CSR(int** _row_id, int** _col, double** _value) {
	free(*_row_id); *_row_id = NULL;
	free(*_col); *_col = NULL;
	free(*_value); *_value = NULL;
}

void COO_to_CSR(int N, int M, int nz, double* val, int* I, int* J, int** _row_id, int** _col, double** _value) {
	int* row_id = (int*)malloc((N + 1) * sizeof(int));
	int* col = (int*)malloc(nz * sizeof(int));
	double* value = (double*)malloc(nz * sizeof(double));
	memset(row_id, 0, sizeof(int) * (N + 1));
	for (int i = 0; i < nz; ++i) {
		row_id[I[i] + 1]++;
	}
	for (int i = 1; i < N + 1; ++i) {
		row_id[i] += row_id[i - 1];
	}
	for (int i = 0; i < nz; ++i) {
		col[row_id[I[i]]] = J[i];
		value[row_id[I[i]]] = val[i];
		++row_id[I[i]];
	}
	for (int i = N; i - 1 >= 0; --i) {
		row_id[i] = row_id[i - 1];
	}
	row_id[0] = 0;
	*_row_id = row_id;
	*_col = col;
	*_value = value;
}
