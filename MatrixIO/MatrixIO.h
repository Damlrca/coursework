#ifndef MATRIX_IO_H
#define MATRIX_IO_H

#include "MatrixUtils.h"

// read unsymmetric or symmetrix mtx matrix from file to COO format
int read_matrix_MTX(const char* fname, matrix_COO* mtx);

// save matrix in COO format to BIN file
int save_matrix_BIN(const char* fname, matrix_COO* mtx);

// read matrix from BIN file to COO format
int read_matrix_BIN(const char* fname, matrix_COO* mtx);

int __read_matrix_MTX(const char* fname, int* M_, int* N_, int* nz_, double** val_, int** I_, int** J_);

int __save_matrix_BIN(const char* fname, int M_, int N_, int nz_, double* val_, int* I_, int* J_);

int __read_matrix_BIN(const char* fname, int* M_, int* N_, int* nz_, double** val_, int** I_, int** J_);

#endif // !MATRIX_IO_H
