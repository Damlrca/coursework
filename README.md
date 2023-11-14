## coursework repo

Умножение разреженных матриц

#### MatrixMult

Несколько функций умножения разреженных матриц

#### MatrixMultTests

Тестирование функций из **MatrixMult**

#### MatrixTests

Тестирование различных функций из **MatrixUtils** и **MatrixIO**

#### MatrixUtils

Стуктуры и функции, необходимые для работы с разреженными матрицами

- matrix_COO // structure for matrix in COO format
- matrix_CSR // structure for matrix in CSR format
- void delete_COO(matrix_COO* mtx) // delete matrix in COO format
- void delete_CSR(matrix_CSR* mtx) // delete matrix in CSR format
- void move_CSR(matrix_CSR* from, matrix_CSR* to) // move one matrix in CSR format from one object to another
- void COO_to_CSR(matrix_COO* mtx_coo, matrix_CSR* mtx_csr) // convert COO to CSR
- void create_transposed(matrix_CSR* from, matrix_CSR* to) // create transposed matrix in CSR format from matrix in CSR format
- void transpose_this(matrix_CSR* mtx) // transpose this matrix in CSR format
- void matrix_addition(matrix_CSR* mtx1, matrix_CSR* mtx2, matrix_CSR* res) // addition of two ORDERED!!! CSR matrices

#### MatrixIO

Функции для чтения и сохранения разреженных матриц в форматах MTX и BIN

- int read_matrix_MTX(const char* fname, matrix_COO* mtx) // read unsymmetric or symmetrix mtx matrix from file to COO format
- int read_matrix_BIN(const char* fname, matrix_COO* mtx) // read matrix from BIN file to COO format
- int save_matrix_BIN(const char* fname, matrix_COO* mtx) // save matrix in COO format to BIN file

#### MatrixGenerator

Генерация матриц для тестов

#### OtherUtils

Другие полезные утилиты

#### sample_matrices

Матрицы для тестов

