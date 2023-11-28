#ifndef MATRIX_MULT_H
#define MATRIX_MULT_H

extern "C" {
#include "MatrixUtils.h"
}

int matrix_mult_naive(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr);

int matrix_mult_naive_2(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr);

#endif // !MATRIX_MULT_H
