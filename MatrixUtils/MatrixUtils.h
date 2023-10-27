#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

// delete matrix int COO format
void delete_COO(double** val, int** I, int** J);

// delete matrix int CSR format
void delete_CSR(int** _row_id, int** _col, double** _value);

// convert COO to CSR
void COO_to_CSR(int N, int M, int nz, double* val, int* I, int* J, int** _row_id, int** _col, double** _value);

// create transposed matrix in CSR format from matrix in CSR format
void create_transposed(int N, int M, int* row_id, int* col, double* value,
	int* N_T, int* M_T, int** row_id_T, int** col_T, double** value_T);

#endif // !MATRIX_UTILS_H
