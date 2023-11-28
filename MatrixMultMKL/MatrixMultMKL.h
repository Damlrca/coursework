#ifndef MATRIX_MULT_MKL_H
#define MATRIX_MULT_MKL_H

#include "MatrixUtils.h"

int matrix_mult_mkl(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* res_csr);

#endif // MATRIX_MULT_MKL_H
