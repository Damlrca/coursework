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

void basic_test(matrix_CSR& matrA_CSR, matrix_CSR& matrB_CSR) {
	// matrix_mult_mkl

	//cout << "C = A * B" << endl;
	matrix_CSR matrC_CSR;
	test_algo("matrix_mult_mkl", matrix_mult_mkl, matrA_CSR, matrB_CSR, matrC_CSR);
	//cout << "Matrix C: N: " << matrC_CSR.N << ", M: " << matrC_CSR.M << ", nz: " << matrC_CSR.row_id[matrC_CSR.N] << endl;

	transpose_this_CSR(&matrC_CSR);
	transpose_this_CSR(&matrC_CSR);
	cout << endl;

	transpose_this_CSR(&matrA_CSR); // !!! matrA_CSR should be sorted
	transpose_this_CSR(&matrA_CSR);
	transpose_this_CSR(&matrB_CSR); // !!! transpose matrix B !!!
	/*
	// matrC_CSR_naive_1
	for (int i = 0; i < 3; i++) {
		matrix_CSR matrC_CSR_naive_1;
		//cout << "C_naive = A * B" << endl;
		test_algo("matrix_mult_naive_1", matrix_mult_naive_1, matrA_CSR, matrB_CSR, matrC_CSR_naive_1);
		//cout << "Matrix C_naive_1: N: " << matrC_CSR_naive_1.N << ", M: " << matrC_CSR_naive_1.M << ", nz: " << matrC_CSR_naive_1.row_id[matrC_CSR_naive_1.N] << endl;

		transpose_this_CSR(&matrC_CSR_naive_1);
		transpose_this_CSR(&matrC_CSR_naive_1);
		comp_algo("matrC_CSR", "matrC_CSR_naive_1", matrC_CSR, matrC_CSR_naive_1);
		delete_CSR(&matrC_CSR_naive_1);
	}
	cout << endl;

	// matrC_CSR_naive_2
	for (int i = 0; i < 3; i++) {
		matrix_CSR matrC_CSR_naive_2;
		//cout << "C_naive_2 = A * B" << endl;
		test_algo("matrix_mult_naive_2", matrix_mult_naive_2, matrA_CSR, matrB_CSR, matrC_CSR_naive_2);
		//cout << "Matrix C_naive_2: N: " << matrC_CSR_naive_2.N << ", M: " << matrC_CSR_naive_2.M << ", nz: " << matrC_CSR_naive_2.row_id[matrC_CSR_naive_2.N] << endl;

		transpose_this_CSR(&matrC_CSR_naive_2);
		transpose_this_CSR(&matrC_CSR_naive_2);
		comp_algo("matrC_CSR", "matrC_CSR_naive_2", matrC_CSR, matrC_CSR_naive_2);
		delete_CSR(&matrC_CSR_naive_2);
	}
	cout << endl;

	// matrC_CSR_naive_3
	for (int i = 0; i < 3; i++) {
		matrix_CSR matrC_CSR_naive_3;
		//cout << "C_naive_3 = A * B" << endl;
		test_algo("matrix_mult_naive_3", matrix_mult_naive_3, matrA_CSR, matrB_CSR, matrC_CSR_naive_3);
		//cout << "Matrix C_naive_3: N: " << matrC_CSR_naive_3.N << ", M: " << matrC_CSR_naive_3.M << ", nz: " << matrC_CSR_naive_3.row_id[matrC_CSR_naive_3.N] << endl;

		transpose_this_CSR(&matrC_CSR_naive_3);
		transpose_this_CSR(&matrC_CSR_naive_3);
		comp_algo("matrC_CSR", "matrC_CSR_naive_3", matrC_CSR, matrC_CSR_naive_3);
		delete_CSR(&matrC_CSR_naive_3);
	}
	cout << endl;
	*/
	// matrC_CSR_naive_1_naiveomp
	for (int i = 0; i < 3; i++) {
		matrix_CSR matrC_CSR_naive_1_naiveomp;
		//cout << "C_naive_1_naiveomp = A * B" << endl;
		test_algo("matrix_mult_naive_1_naiveomp", matrix_mult_naive_1_naiveomp, matrA_CSR, matrB_CSR, matrC_CSR_naive_1_naiveomp);
		//cout << "Matrix C_naive_1_naiveomp: N: " << matrC_CSR_naive_1_naiveomp.N << ", M: " << matrC_CSR_naive_1_naiveomp.M << ", nz: " << matrC_CSR_naive_1_naiveomp.row_id[matrC_CSR_naive_1_naiveomp.N] << endl;

		transpose_this_CSR(&matrC_CSR_naive_1_naiveomp);
		transpose_this_CSR(&matrC_CSR_naive_1_naiveomp);
		comp_algo("matrC_CSR", "matrC_CSR_naive_1_naiveomp", matrC_CSR, matrC_CSR_naive_1_naiveomp);
		delete_CSR(&matrC_CSR_naive_1_naiveomp);
	}
	cout << endl;

	// matrC_CSR_naive_2_naiveomp
	for (int i = 0; i < 3; i++) {
		matrix_CSR matrC_CSR_naive_2_naiveomp;
		//cout << "C_naive_2_naiveomp = A * B" << endl;
		test_algo("matrix_mult_naive_2_naiveomp", matrix_mult_naive_2_naiveomp, matrA_CSR, matrB_CSR, matrC_CSR_naive_2_naiveomp);
		//cout << "Matrix C_naive_2_naiveomp: N: " << matrC_CSR_naive_2_naiveomp.N << ", M: " << matrC_CSR_naive_2_naiveomp.M << ", nz: " << matrC_CSR_naive_2_naiveomp.row_id[matrC_CSR_naive_2_naiveomp.N] << endl;

		transpose_this_CSR(&matrC_CSR_naive_2_naiveomp);
		transpose_this_CSR(&matrC_CSR_naive_2_naiveomp);
		comp_algo("matrC_CSR", "matrC_CSR_naive_2_naiveomp", matrC_CSR, matrC_CSR_naive_2_naiveomp);
		delete_CSR(&matrC_CSR_naive_2_naiveomp);
	}
	cout << endl;

	// matrC_CSR_naive_3_naiveomp
	for (int i = 0; i < 3; i++) {
		matrix_CSR matrC_CSR_naive_3_naiveomp;
		//cout << "C_naive_3_naiveomp = A * B" << endl;
		test_algo("matrix_mult_naive_3_naiveomp", matrix_mult_naive_3_naiveomp, matrA_CSR, matrB_CSR, matrC_CSR_naive_3_naiveomp);
		//cout << "Matrix C_naive_3_naiveomp: N: " << matrC_CSR_naive_3_naiveomp.N << ", M: " << matrC_CSR_naive_3_naiveomp.M << ", nz: " << matrC_CSR_naive_3_naiveomp.row_id[matrC_CSR_naive_3_naiveomp.N] << endl;

		transpose_this_CSR(&matrC_CSR_naive_3_naiveomp);
		transpose_this_CSR(&matrC_CSR_naive_3_naiveomp);
		comp_algo("matrC_CSR", "matrC_CSR_naive_3_naiveomp", matrC_CSR, matrC_CSR_naive_3_naiveomp);
		delete_CSR(&matrC_CSR_naive_3_naiveomp);
	}
	cout << endl;

	// matrC_CSR_naive_1_queueomp
	for (int i = 0; i < 3; i++) {
		matrix_CSR matrC_CSR_naive_1_queueomp;
		//cout << "C_naive_1_queueomp = A * B" << endl;
		test_algo("matrix_mult_naive_1_queueomp", matrix_mult_naive_1_queueomp, matrA_CSR, matrB_CSR, matrC_CSR_naive_1_queueomp);
		//cout << "Matrix C_naive_1_queueomp: N: " << matrC_CSR_naive_1_queueomp.N << ", M: " << matrC_CSR_naive_1_queueomp.M << ", nz: " << matrC_CSR_naive_1_queueomp.row_id[matrC_CSR_naive_1_queueomp.N] << endl;

		transpose_this_CSR(&matrC_CSR_naive_1_queueomp);
		transpose_this_CSR(&matrC_CSR_naive_1_queueomp);
		comp_algo("matrC_CSR", "matrC_CSR_naive_1_queueomp", matrC_CSR, matrC_CSR_naive_1_queueomp);
		delete_CSR(&matrC_CSR_naive_1_queueomp);
	}
	cout << endl;

	// matrC_CSR_naive_2_queueomp
	for (int i = 0; i < 3; i++) {
		matrix_CSR matrC_CSR_naive_2_queueomp;
		//cout << "C_naive_2_queueomp = A * B" << endl;
		test_algo("matrix_mult_naive_2_queueomp", matrix_mult_naive_2_queueomp, matrA_CSR, matrB_CSR, matrC_CSR_naive_2_queueomp);
		//cout << "Matrix C_naive_2_queueomp: N: " << matrC_CSR_naive_2_queueomp.N << ", M: " << matrC_CSR_naive_2_queueomp.M << ", nz: " << matrC_CSR_naive_2_queueomp.row_id[matrC_CSR_naive_2_queueomp.N] << endl;

		transpose_this_CSR(&matrC_CSR_naive_2_queueomp);
		transpose_this_CSR(&matrC_CSR_naive_2_queueomp);
		comp_algo("matrC_CSR", "matrC_CSR_naive_2_queueomp", matrC_CSR, matrC_CSR_naive_2_queueomp);
		delete_CSR(&matrC_CSR_naive_2_queueomp);
	}
	cout << endl;

	// matrC_CSR_naive_3_queueomp
	for (int i = 0; i < 3; i++) {
		matrix_CSR matrC_CSR_naive_3_queueomp;
		//cout << "C_naive_3_queueomp = A * B" << endl;
		test_algo("matrix_mult_naive_3_queueomp", matrix_mult_naive_3_queueomp, matrA_CSR, matrB_CSR, matrC_CSR_naive_3_queueomp);
		//cout << "Matrix C_naive_3_queueomp: N: " << matrC_CSR_naive_3_queueomp.N << ", M: " << matrC_CSR_naive_3_queueomp.M << ", nz: " << matrC_CSR_naive_3_queueomp.row_id[matrC_CSR_naive_3_queueomp.N] << endl;

		transpose_this_CSR(&matrC_CSR_naive_3_queueomp);
		transpose_this_CSR(&matrC_CSR_naive_3_queueomp);
		comp_algo("matrC_CSR", "matrC_CSR_naive_3_queueomp", matrC_CSR, matrC_CSR_naive_3_queueomp);
		delete_CSR(&matrC_CSR_naive_3_queueomp);
	}
	cout << endl;

	delete_CSR(&matrC_CSR);
}

int main() {
	cout << "-------------------------------------------" << endl;
	cout << "test uniform matrix 10'000 x 10'000, k = 15" << endl;
	cout << "-------------------------------------------" << endl;
	cout << endl;

	matrix_CSR matrA_CSR;
	generate_uniform_square_sparse_matrix_CSR(10'000, 15, matrA_CSR);
	cout << "Matrix A generated successfully" << endl;
	cout << "Matrix A: N: " << matrA_CSR.N << ", M: " << matrA_CSR.M << ", nz: " << matrA_CSR.row_id[matrA_CSR.N] << endl;
	matrix_CSR matrB_CSR;
	generate_uniform_square_sparse_matrix_CSR(10'000, 15, matrB_CSR);
	cout << "Matrix B generated successfully" << endl;
	cout << "Matrix B: N: " << matrB_CSR.N << ", M: " << matrB_CSR.M << ", nz: " << matrB_CSR.row_id[matrB_CSR.N] << endl;
	cout << endl;

	basic_test(matrA_CSR, matrB_CSR);

	delete_CSR(&matrA_CSR);
	//delete_CSR(&matrB_CSR);

	cout << "-------------------------------------------" << endl;
	cout << "test nonuniform matrix 10'000 x 10'000, k1 = 30, k2 = 1" << endl;
	cout << "-------------------------------------------" << endl;
	cout << endl;

	generate_nonuniform_square_sparse_matrix_CSR(10'000, 30, 1, matrA_CSR);
	cout << "Matrix A generated successfully" << endl;
	cout << "Matrix A: N: " << matrA_CSR.N << ", M: " << matrA_CSR.M << ", nz: " << matrA_CSR.row_id[matrA_CSR.N] << endl;
	//generate_nonuniform_square_sparse_matrix_CSR(10'000, 30, 1, matrB_CSR);
	//cout << "Matrix B generated successfully" << endl;
	//cout << "Matrix B: N: " << matrB_CSR.N << ", M: " << matrB_CSR.M << ", nz: " << matrB_CSR.row_id[matrB_CSR.N] << endl;
	transpose_this_CSR(&matrB_CSR);
	cout << endl;

	basic_test(matrA_CSR, matrB_CSR);

	delete_CSR(&matrA_CSR);
	delete_CSR(&matrB_CSR);

	return 0;
}
