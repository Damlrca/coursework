#ifndef MATRIX_GENERATOR_HPP
#define MATRIX_GENERATOR_HPP

extern "C" {
#include "MatrixUtils.h"
}

// matrix "res" should be empty!
// ��������� ���������� ������� �������� n*n �� k �������� � ������ ������
void generate_uniform_square_sparse_matrix_CSR(int size, int k, matrix_CSR& res);

// matrix "res" should be empty!
// ��������� ���������� ������� �������� n*n
// ���-�� ��������� � ������� ������� �� k1 �� k2
void generate_nonuniform_square_sparse_matrix_CSR(int size, int k1, int k2, matrix_CSR& res);

#endif // !MATRIX_GENERATOR_HPP
