#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

typedef struct _matrix_COO {
	int N, M, nz;
	double* val;
	int* I;
	int* J;
} matrix_COO;

typedef struct _matrix_CSR {
	int N, M;
	int* row_id;
	int* col;
	double* value;
} matrix_CSR;

// delete matrix in COO format
void delete_COO(matrix_COO* mtx);

// delete matrix in CSR format
void delete_CSR(matrix_CSR* mtx);

// convert COO to CSR
void COO_to_CSR(matrix_COO* mtx_coo, matrix_CSR* mtx_csr);
void _COO_to_CSR(int N, int M, int nz, double* val, int* I, int* J, int** _row_id, int** _col, double** _value);

// create transposed matrix in CSR format from matrix in CSR format
void create_transposed(int N, int M, int* row_id, int* col, double* value,
	int* N_T, int* M_T, int** row_id_T, int** col_T, double** value_T);

// transpose this matrix in CSR format
void transpose_this(int* N, int* M, int** row_id, int** col, double** value);

// addition of two ORDERED!!! CSR matrices
void matrix_addition(int N1, int M1, int* row_id1, int* col1, double* value1,
	int N2, int M2, int* row_id2, int* col2, double* value2,
	int* N, int* M, int** row_id, int** col, double** value);

#endif // !MATRIX_UTILS_H
