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

// move one matrix in CSR format from one object to another
void move_CSR(matrix_CSR* from, matrix_CSR* to);

// convert COO to CSR
void COO_to_CSR(matrix_COO* mtx_coo, matrix_CSR* mtx_csr);
void _COO_to_CSR(int N, int M, int nz, double* val, int* I, int* J, int** _row_id, int** _col, double** _value);

// create transposed matrix in CSR format from matrix in CSR format
void create_transposed(matrix_CSR* from, matrix_CSR* to);
void _create_transposed(int N, int M, int* row_id, int* col, double* value,
	int* N_T, int* M_T, int** row_id_T, int** col_T, double** value_T);

// transpose this matrix in CSR format
void transpose_this(matrix_CSR* mtx);

// addition of two ORDERED!!! CSR matrices
void matrix_addition(matrix_CSR* mtx1, matrix_CSR* mtx2, matrix_CSR* res);

#endif // !MATRIX_UTILS_H
