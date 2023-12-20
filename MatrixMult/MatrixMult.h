#ifndef MATRIX_MULT_H
#define MATRIX_MULT_H

extern "C" {
#include "MatrixUtils.h"
}

// Пока что только квадратные матрицы!!!

// Матрица B должна быть транспонирована!
// Первый, самый медленный алгоритм, не используем то, что элементы в строках отсортированы по индексам столбцов
int matrix_mult_naive_1(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr);

// Матрица B должна быть транспонирована!
// В матрицах A и B элементы в строках должны быть отсортированы по индексам стобцов
int matrix_mult_naive_2(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr);

// Матрица B должна быть транспонирована!
// "остроумный вариант вычисления скалярного произведения разреженных веторов"
int matrix_mult_naive_3(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr);

int matrix_mult_naive_1_naiveomp(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr);
int matrix_mult_naive_2_naiveomp(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr);
int matrix_mult_naive_3_naiveomp(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr);

int matrix_mult_naive_1_queueomp(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr);
int matrix_mult_naive_2_queueomp(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr);
int matrix_mult_naive_3_queueomp(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr);

#endif // !MATRIX_MULT_H
