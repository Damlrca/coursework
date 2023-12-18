## coursework repo

Умножение разреженных матриц!

#### MatrixMult (C, library)

Несколько функций умножения разреженных матриц

#### MatrixMultMKL (C, library)

Функция умножения разреженных матриц с использованием функции из oneMKL

- int **matrix_mult_mkl**(matrix_CSR\* A_csr, matrix_CSR\* B_csr, matrix_CSR\* res_csr);

#### MatrixMultTests (Cpp, executable)

Тестирование функций из **MatrixMult**, **MatrixMultMKL**

#### MatrixTests (Cpp, executable)

Тестирование различных функций из **MatrixUtils**, **MatrixIO**, **OtherUtils**

#### MatrixUtils (C, library)

Структуры и функции, необходимые для работы с разреженными матрицами

- **matrix_COO** : structure for sparse matrix in COO format
- **matrix_CSR** : structure for sparse matrix in CSR format
- void **delete_COO**(matrix_COO\* mtx) : delete contents of COO-matrix
- void **delete_CSR**(matrix_CSR\* mtx) : delete contents of CSR-matrix
- void **move_CSR**(matrix_CSR\* from, matrix_CSR\* to) : move contents of one CSR-matrix to another
- void **convert_COO_to_CSR**(matrix_COO\* mtx_coo, matrix_CSR\* mtx_csr) : convert COO-matrix into CSR-matrix
- void **create_transposed_CSR**(matrix_CSR\* from, matrix_CSR\* to) : create transposed CSR-matrix from CSR-matrix
- void **transpose_this_CSR**(matrix_CSR\* mtx) : transpose this CSR-matrix
- void **matrix_addition**(matrix_CSR\* mtx1, matrix_CSR\* mtx2, matrix_CSR\* res) : addition of two CSR-matrices

#### MatrixIO (C, library)

Функции для чтения и сохранения разреженных матриц в форматах MTX и BIN

- int **read_matrix_MTX**(const char\* fname, matrix_COO\* mtx) : read unsymmetric or symmetrix matrix from mtx file to COO-matrix
- int **read_matrix_BIN**(const char\* fname, matrix_COO\* mtx) : read matrix from BIN file to COO-matrix
- int **save_matrix_BIN**(const char\* fname, matrix_COO\* mtx) : save COO-matrix to BIN file

#### MatrixGenerator (Cpp, library)

Генерация матриц для тестов

- void **generate_uniform_square_sparse_matrix_CSR**(int size, int k, matrix_CSR& res) : генерация квадратной матрицы размером n\*n по k элеметов в каждой строке
- void **generate_nonuniform_square_sparse_matrix_CSR**(int size, int k1, int k2, matrix_CSR& res) : генерация квадратной матрицы размером n\*n кол-во элементов в строках линейно от k1 до k2

#### OtherUtils (Cpp, library)

Другие полезные классы и функции

- **MyTimer** : class for measuring time
	- static void **MyTimer::SetStartTime**() : sets start time of measurment
	- static void **MyTimer::SetEndTime**() : sets end time of measurment
	- static long long **GetDifferenceMs**() : get difference between start and end time in milliseconds
- void **print_dense**(std::vector<std::vector<double>>& dense) : print dense matrix
- void **print_COO_as_dense**(const matrix_COO& mtx, const char\* mtx_name) : print COO-matrix as dense matrix
- void **print_CSR**(const matrix_CSR& mtx) : print CSR-matrix in CSR format
- void **print_CSR_as_dense**(const matrix_CSR& mtx, const char\* mtx_name) : print CSR-matrix as dense matrix
- void **matrix_difference**(matrix_CSR\* mtx1, matrix_CSR\* mtx2, matrix_CSR\* res) : difference of two CSR-matrices
- bool **matrix_compare**(matrix_CSR\* mtx1, matrix_CSR\* mtx2) : returns true if matrices "mtx1" and "mtx2" are equal

#### MTXtoBINconverter (Cpp, executable)

Создание файла матрицы в BIN формате на основе MTX-файла

#### sample_matrices

Матрицы для тестов

