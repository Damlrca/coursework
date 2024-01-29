#ifndef MATRIX_MULT_TWO_PHASES_H
#define MATRIX_MULT_TWO_PHASES_H

extern "C" {
#include "MatrixUtils.h"
}

// Пока что только квадратные матрицы!!!

int matrix_mult_TF_first_phase(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr, int delta);
int matrix_mult_TF_second_phase(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr, int delta);

#endif // !MATRIX_MULT_TWO_PHASES_H
