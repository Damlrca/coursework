#include "MatrixUtils.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void delete_COO(matrix_COO* mtx) {
	mtx->N = 0;
	mtx->M = 0;
	mtx->nz = 0;
	free(mtx->val); mtx->val = NULL;
	free(mtx->I); mtx->I = NULL;
	free(mtx->J); mtx->J = NULL;
}

void delete_CSR(matrix_CSR* mtx) {
	mtx->N = 0;
	mtx->M = 0;
	free(mtx->value); mtx->value = NULL;
	free(mtx->row_id); mtx->row_id = NULL;
	free(mtx->col); mtx->col = NULL;
}

void move_CSR(matrix_CSR* from, matrix_CSR* to) {
	to->N = from->N; from->N = 0;
	to->M = from->M; from->M = 0;
	to->row_id = from->row_id; from->row_id = NULL;
	to->col = from->col; from->col = NULL;
	to->value = from->value; from->value = NULL;
}

void convert_COO_to_CSR(matrix_COO* mtx_coo, matrix_CSR* mtx_csr) {
	_convert_COO_to_CSR(mtx_coo->N, mtx_coo->M, mtx_coo->nz, mtx_coo->val, mtx_coo->I, mtx_coo->J,
		&mtx_csr->row_id, &mtx_csr->col, &mtx_csr->value);
	mtx_csr->N = mtx_coo->N;
	mtx_csr->M = mtx_coo->M;
}

void _convert_COO_to_CSR(int N, int M, int nz, double* val, int* I, int* J, int** _row_id, int** _col, double** _value) {
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

void create_transposed_CSR(matrix_CSR* from, matrix_CSR* to) {
	_create_transposed_CSR(from->N, from->M, from->row_id, from->col, from->value,
		&to->N, &to->M, &to->row_id, &to->col, &to->value);
}

void _create_transposed_CSR(int N, int M, int* row_id, int* col, double* value,
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

void transpose_this_CSR(matrix_CSR* mtx) {
	matrix_CSR res;
	create_transposed_CSR(mtx, &res);
	delete_CSR(mtx);
	move_CSR(&res, mtx);
}

void matrix_addition(matrix_CSR* mtx1, matrix_CSR* mtx2, matrix_CSR* res) {
	int nz = 0;
	int N = mtx1->N > mtx2->N ? mtx1->N : mtx2->N;
	int M = mtx1->M > mtx2->M ? mtx1->M : mtx2->M;
	int* row_id = (int*)malloc((N + 1) * sizeof(int));
	int* col;
	double* value;
	int i = 0;
	int pr = 0;
	memset(row_id, 0, (N + 1) * sizeof(int));
	for (; i < mtx1->N && i < mtx2->N; ++i) {
		int a1 = mtx1->row_id[i], b1 = mtx1->row_id[i + 1];
		int a2 = mtx2->row_id[i], b2 = mtx2->row_id[i + 1];
		while (a1 != b1 && a2 != b2) {
			if (mtx1->col[a1] < mtx2->col[a2]) {
				++row_id[i + 1]; ++a1;
			}
			else if (mtx1->col[a1] > mtx2->col[a2]) {
				++row_id[i + 1]; ++a2;
			}
			else { // mtx1->col[a1] == mtx2->col[a2]
				++row_id[i + 1]; ++a1; ++a2;
			}
		}
		row_id[i + 1] += b1 - a1;
		row_id[i + 1] += b2 - a2;
	}
	for (; i < mtx1->N; ++i) {
		int a1 = mtx1->row_id[i], b1 = mtx1->row_id[i + 1];
		row_id[i + 1] += b1 - a1;
	}
	for (; i < mtx2->N; ++i) {
		int a2 = mtx2->row_id[i], b2 = mtx2->row_id[i + 1];
		row_id[i + 1] += b2 - a2;
	}
	for (int i = 1; i < N + 1; ++i) {
		nz += pr;
		pr = row_id[i];
		row_id[i] = nz;
	}
	nz += pr;
	col = (int*)malloc(nz * sizeof(int));
	value = (double*)malloc(nz * sizeof(double));
	i = 0;
	for (; i < mtx1->N && i < mtx2->N; ++i) {
		int a1 = mtx1->row_id[i], b1 = mtx1->row_id[i + 1];
		int a2 = mtx2->row_id[i], b2 = mtx2->row_id[i + 1];
		while (a1 != b1 && a2 != b2) {
			int Index = row_id[i + 1];
			if (mtx1->col[a1] < mtx2->col[a2]) {
				col[Index] = mtx1->col[a1];
				value[Index] = mtx1->value[a1];
				++a1;
			}
			else if (mtx1->col[a1] > mtx2->col[a2]) {
				col[Index] = mtx2->col[a2];
				value[Index] = mtx2->value[a2];
				++a2;
			}
			else { // mtx1->col[a1] == mtx2->col[a2]
				col[Index] = mtx1->col[a1];
				value[Index] = mtx1->value[a1] + mtx2->value[a2];
				++a1; ++a2;
			}
			++row_id[i + 1];
		}
		while (a1 != b1) {
			int Index = row_id[i + 1];
			col[Index] = mtx1->col[a1];
			value[Index] = mtx1->value[a1];
			++a1;
			++row_id[i + 1];
		}
		while (a2 != b2) {
			int Index = row_id[i + 1];
			col[Index] = mtx2->col[a2];
			value[Index] = mtx2->value[a2];
			++a2;
			++row_id[i + 1];
		}
	}
	for (; i < mtx1->N; ++i) {
		int a1 = mtx1->row_id[i], b1 = mtx1->row_id[i + 1];
		while (a1 != b1) {
			int Index = row_id[i + 1];
			col[Index] = mtx1->col[a1];
			value[Index] = mtx1->value[a1];
			++a1;
			++row_id[i + 1];
		}
	}
	for (; i < mtx2->N; ++i) {
		int a2 = mtx2->row_id[i], b2 = mtx2->row_id[i + 1];
		while (a2 != b2) {
			int Index = row_id[i + 1];
			col[Index] = mtx2->col[a2];
			value[Index] = mtx2->value[a2];
			++a2;
			++row_id[i + 1];
		}
	}
	res->N = N;
	res->M = M;
	res->row_id = row_id;
	res->col = col;
	res->value = value;
}
