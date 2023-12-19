#ifndef MATRIX_MULT_H
#define MATRIX_MULT_H

extern "C" {
#include "MatrixUtils.h"
}

// ���� ��� ������ ���������� �������!!!

// ������� B ������ ���� ���������������!
// ������, ����� ��������� ��������, �� ���������� ��, ��� �������� � ������� ������������� �� �������� ��������
int matrix_mult_naive_1(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr);

// ������� B ������ ���� ���������������!
// � �������� A � B �������� � ������� ������ ���� ������������� �� �������� �������
int matrix_mult_naive_2(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr);

// ������� B ������ ���� ���������������!
// "���������� ������� ���������� ���������� ������������ ����������� �������"
int matrix_mult_naive_3(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr);

int matrix_mult_naive_1_naiveomp(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr);
int matrix_mult_naive_2_naiveomp(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr);
int matrix_mult_naive_3_naiveomp(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr);

int matrix_mult_naive_1_queueomp(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr);
int matrix_mult_naive_2_queueomp(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr);
int matrix_mult_naive_3_queueomp(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr);

#endif // !MATRIX_MULT_H
