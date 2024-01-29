#ifndef MATRIX_GENERATOR_HPP
#define MATRIX_GENERATOR_HPP

extern "C" {
#include "MatrixUtils.h"
}

// matrix "res" should be empty!
// генерация квадратной матрицы размером n*n по k элеметов в каждой строке
void generate_uniform_square_sparse_matrix(int n, int k, matrix_CSR& res);

// matrix "res" should be empty!
// генерация квадратной матрицы размером n*n
// кол-во элементов в строках линейно от k1 до k2
void generate_nonuniform_square_sparse_matrix(int n, int k1, int k2, matrix_CSR& res);

#endif // !MATRIX_GENERATOR_HPP
