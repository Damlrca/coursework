#ifndef MATRIX_MULT_H
#define MATRIX_MULT_H

extern "C" {
#include "MatrixUtils.h"
}

// ���� ��� ������ ���������� �������!!!

// ������� B ������ ���� ���������������!
// ������, ����� ��������� ��������, �� ���������� ��, ��� �������� � ������� ������������� �� �������� ��������
int matrix_mult_naive(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr);

// ������� B ������ ���� ���������������!
// � �������� A � B �������� � ������� ������ ���� ������������� �� �������� �������
int matrix_mult_naive_2(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr);

// ������� B ������ ���� ���������������!
// "���������� ������� ���������� ���������� ������������ ����������� �������"
int matrix_mult_naive_3(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr);

#endif // !MATRIX_MULT_H
