#include "MatrixUtils.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void delete_COO(double** val, int** I, int** J) {
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
	int S = 0, pr = 0;
	memset(row_id, 0, (N + 1) * sizeof(int));
	for (int i = 0; i < nz; ++i) {
		++row_id[I[i] + 1];
	}
	for (int i = 1; i < N + 1; ++i) {
		S += pr;
		pr = row_id[i];
		row_id[i] = S;
	}
	for (int i = 0; i < nz; ++i) {
		int RIndex = row_id[I[i] + 1];
		col[RIndex] = J[i];
		value[RIndex] = val[i];
		++row_id[I[i] + 1];
	}
	*_row_id = row_id;
	*_col = col;
	*_value = value;
}

void create_transposed(int N, int M, int* row_id, int* col, double* value,
	int* N_T, int* M_T, int** _row_id_T, int** _col_T, double** _value_T) {
	int nz = row_id[N];
	*N_T = M;
	*M_T = N;
	int* row_id_T = (int*)malloc((*N_T + 1) * sizeof(int));
	int* col_T = (int*)malloc(nz * sizeof(int));
	double* value_T = (double*)malloc(nz * sizeof(double));
	int S = 0, pr = 0;
	memset(row_id_T, 0, (*N_T + 1) * sizeof(int));
	for (int i = 0; i < nz; ++i) {
		++row_id_T[col[i] + 1];
	}
	for (int i = 1; i < *N_T + 1; ++i) {
		S += pr;
		pr = row_id_T[i];
		row_id_T[i] = S;
	}
	for (int i = 0; i < N; i++) {
		int a = row_id[i];
		int b = row_id[i + 1];
		while (a != b) {
			int RIndex = row_id_T[col[a] + 1];
			col_T[RIndex] = i;
			value_T[RIndex] = value[a];
			++row_id_T[col[a] + 1];
			++a;
		}
	}
	*_row_id_T = row_id_T;
	*_col_T = col_T;
	*_value_T = value_T;
}
