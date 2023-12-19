extern "C" {
#include "MatrixUtils.h"
#include "MatrixIO.h"
#include "MatrixMultMKL.h"
}
#include "MatrixGenerator.hpp"
#include "MatrixMult.h"
#include "OtherUtils.hpp"
#include <iostream>

using std::cout;
using std::endl;

//const char* fileA = "D:\\source\\coursework\\sample_matrices\\multA.bin";
//const char* fileB = "D:\\source\\coursework\\sample_matrices\\multB.bin";

void test_algo(const char* algo_name, int (*algo) (matrix_CSR*, matrix_CSR*, matrix_CSR*), matrix_CSR& A_csr, matrix_CSR& B_csr, matrix_CSR& C_csr) {
	MyTimer::SetStartTime();
	int status = algo(&A_csr, &B_csr, &C_csr);
	MyTimer::SetEndTime();
	if (status != 0) {
		cout << "(ERROR) : " << algo_name << " failed" << endl;
		exit(-1);
	}
	cout << "(OK) : " << "matrix multiplied using \"" << algo_name << "\" succesfully in " << MyTimer::GetDifferenceMs() << "ms" << endl;
}

void comp_algo(const char* algo_name_1, const char* algo_name_2, matrix_CSR& C_csr_1, matrix_CSR& C_csr_2) {
	if (matrix_compare(&C_csr_1, &C_csr_2)) {
		cout << "(OK) : " << algo_name_1 << " == " << algo_name_2 << endl;
	}
	else {
		cout << "(ERROR) : " << algo_name_1 << " != " << algo_name_2 << endl;
		exit(-1);
	}
}

int main() {
	//matrix_COO matrA_COO;
	//read_matrix_BIN(fileA, &matrA_COO);
	//matrix_COO matrB_COO;
	//read_matrix_BIN(fileB, &matrB_COO);
	matrix_CSR matrA_CSR;
	//convert_COO_to_CSR(&matrA_COO, &matrA_CSR);
	//generate_uniform_square_sparse_matrix_CSR(20'000, 25, matrA_CSR);
	generate_uniform_square_sparse_matrix_CSR(10'000, 15, matrA_CSR);
	//generate_nonuniform_square_sparse_matrix_CSR(20'000, 0, 50, matrA_CSR);
	//cout << "A (" << fileA << ")" << endl;
	cout << "Matrix A generated successfully" << endl;
	cout << "Matrix A: N: " << matrA_CSR.N << ", M: " << matrA_CSR.M << ", nz: " << matrA_CSR.row_id[matrA_CSR.N] << endl;
	matrix_CSR matrB_CSR;
	//convert_COO_to_CSR(&matrB_COO, &matrB_CSR);
	//generate_uniform_square_sparse_matrix_CSR(20'000, 25, matrB_CSR);
	generate_uniform_square_sparse_matrix_CSR(10'000, 15, matrB_CSR);
	//generate_nonuniform_square_sparse_matrix_CSR(20'000, 0, 50, matrB_CSR);
	//cout << "B (" << fileB << ")" << endl;
	cout << "Matrix B generated successfully" << endl;
	cout << "Matrix B: N: " << matrB_CSR.N << ", M: " << matrB_CSR.M << ", nz: " << matrB_CSR.row_id[matrB_CSR.N] << endl;
	cout << endl;

	// matrix_mult_mkl

	cout << "C = A * B" << endl;
	matrix_CSR matrC_CSR;
	test_algo("matrix_mult_mkl", matrix_mult_mkl, matrA_CSR, matrB_CSR, matrC_CSR);
	cout << "Matrix C: N: " << matrC_CSR.N << ", M: " << matrC_CSR.M << ", nz: " << matrC_CSR.row_id[matrC_CSR.N] << endl;

	transpose_this_CSR(&matrC_CSR);
	transpose_this_CSR(&matrC_CSR);
	cout << endl;

	//

	transpose_this_CSR(&matrB_CSR); // !!! transpose matrix B !!!

	// matrC_CSR_naive_3

	matrix_CSR matrC_CSR_naive_3;
	cout << "C_naive_3 = A * B" << endl;
	test_algo("matrix_mult_naive_3", matrix_mult_naive_3, matrA_CSR, matrB_CSR, matrC_CSR_naive_3);
	cout << "Matrix C_naive_3: N: " << matrC_CSR_naive_3.N << ", M: " << matrC_CSR_naive_3.M << ", nz: " << matrC_CSR_naive_3.row_id[matrC_CSR_naive_3.N] << endl;

	transpose_this_CSR(&matrC_CSR_naive_3);
	transpose_this_CSR(&matrC_CSR_naive_3);
	comp_algo("matrC_CSR", "matrC_CSR_naive_3", matrC_CSR, matrC_CSR_naive_3);
	cout << endl;

	// matrC_CSR_naive

	matrix_CSR matrC_CSR_naive;
	cout << "C_naive = A * B" << endl;
	test_algo("matrix_mult_naive", matrix_mult_naive, matrA_CSR, matrB_CSR, matrC_CSR_naive);
	cout << "Matrix C_naive: N: " << matrC_CSR_naive.N << ", M: " << matrC_CSR_naive.M << ", nz: " << matrC_CSR_naive.row_id[matrC_CSR_naive.N] << endl;
	
	transpose_this_CSR(&matrC_CSR_naive);
	transpose_this_CSR(&matrC_CSR_naive);
	comp_algo("matrC_CSR", "matrC_CSR_naive", matrC_CSR, matrC_CSR_naive);
	cout << endl;

	transpose_this_CSR(&matrA_CSR); // !!! matrA_CSR should be sorted
	transpose_this_CSR(&matrA_CSR);
	// !!! (matrB_CSR is already transposed and sorted)

	// matrC_CSR_naive_2

	matrix_CSR matrC_CSR_naive_2;
	cout << "C_naive_2 = A * B" << endl;
	test_algo("matrix_mult_naive_2", matrix_mult_naive_2, matrA_CSR, matrB_CSR, matrC_CSR_naive_2);
	cout << "Matrix C_naive_2: N: " << matrC_CSR_naive_2.N << ", M: " << matrC_CSR_naive_2.M << ", nz: " << matrC_CSR_naive_2.row_id[matrC_CSR_naive_2.N] << endl;

	transpose_this_CSR(&matrC_CSR_naive_2);
	transpose_this_CSR(&matrC_CSR_naive_2);
	comp_algo("matrC_CSR", "matrC_CSR_naive_2", matrC_CSR, matrC_CSR_naive_2);
	cout << endl;

	//

	delete_CSR(&matrA_CSR);
	delete_CSR(&matrB_CSR);
	delete_CSR(&matrC_CSR);
	delete_CSR(&matrC_CSR_naive);
	delete_CSR(&matrC_CSR_naive_2);
	delete_CSR(&matrC_CSR_naive_3);

	return 0;
}
