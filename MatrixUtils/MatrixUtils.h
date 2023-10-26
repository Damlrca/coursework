#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

// delete matrix int COO format
void delete_matrix(int* N, int* M, int* nz, double** val, int** I, int** J);

// delete matrix int CSR format
void delete_CSR(int** _row_id, int** _col, double** _value);

// convert COO to CSR
void COO_to_CSR(int N, int M, int nz, double* val, int* I, int* J, int** _row_id, int** _col, double** _value);

#endif // !MATRIX_UTILS_H
