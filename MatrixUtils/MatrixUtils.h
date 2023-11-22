#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

// structure for sparse matrix in COO format
typedef struct _matrix_COO {
	int N, M, nz;
	double* val;
	int* I;
	int* J;
} matrix_COO;

// structure for sparse matrix in CSR format
typedef struct _matrix_CSR {
	int N, M;
	int* row_id;
	int* col;
	double* value;
} matrix_CSR;

// delete contents of COO-matrix
void delete_COO(matrix_COO* mtx);

// delete contents of CSR-matrix
void delete_CSR(matrix_CSR* mtx);

// move contents of one CSR-matrix to another
// matrix "to" should be empty!
void move_CSR(matrix_CSR* from, matrix_CSR* to);

// convert COO-matrix into CSR-matrix
// matrix "mtx_csr" should be empty!
void convert_COO_to_CSR(matrix_COO* mtx_coo, matrix_CSR* mtx_csr);
void _convert_COO_to_CSR(int N, int M, int nz, double* val, int* I, int* J, int** _row_id, int** _col, double** _value);

// create transposed CSR-matrix from CSR-matrix
// matrix "to" should be empty!
void create_transposed_CSR(matrix_CSR* from, matrix_CSR* to);
void _create_transposed_CSR(int N, int M, int* row_id, int* col, double* value,
	int* N_T, int* M_T, int** row_id_T, int** col_T, double** value_T);

// transpose this CSR-matrix
void transpose_this_CSR(matrix_CSR* mtx);

// addition of two CSR-matrices
// matrices "mtx1" and "mtx2" should be ordered!
// matrix "res" should be empty!
void matrix_addition(matrix_CSR* mtx1, matrix_CSR* mtx2, matrix_CSR* res);

#endif // !MATRIX_UTILS_H
