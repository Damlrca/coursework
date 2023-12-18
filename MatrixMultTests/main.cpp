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

int main() {
	//matrix_COO matrA_COO;
	//read_matrix_BIN(fileA, &matrA_COO);
	//matrix_COO matrB_COO;
	//read_matrix_BIN(fileB, &matrB_COO);
	matrix_CSR matrA_CSR;
	//convert_COO_to_CSR(&matrA_COO, &matrA_CSR);
	//generate_uniform_square_sparse_matrix_CSR(20'000, 25, matrA_CSR);
	generate_nonuniform_square_sparse_matrix_CSR(20'000, 0, 50, matrA_CSR);
	//cout << "A (" << fileA << ")" << endl;
	cout << "Matrix A:" << endl;
	cout << "N: " << matrA_CSR.N << ", M: " << matrA_CSR.M << ", nz: " << matrA_CSR.row_id[matrA_CSR.N] << endl;
	matrix_CSR matrB_CSR;
	//convert_COO_to_CSR(&matrB_COO, &matrB_CSR);
	//generate_uniform_square_sparse_matrix_CSR(20'000, 25, matrB_CSR);
	generate_nonuniform_square_sparse_matrix_CSR(20'000, 0, 50, matrB_CSR);
	//cout << "B (" << fileB << ")" << endl;
	cout << "Matrix B:" << endl;
	cout << "N: " << matrB_CSR.N << ", M: " << matrB_CSR.M << ", nz: " << matrB_CSR.row_id[matrB_CSR.N] << endl;

	cout << "C = A * B" << endl;
	matrix_CSR matrC_CSR;
	MyTimer::SetStartTime();
	int status = matrix_mult_mkl(&matrA_CSR, &matrB_CSR, &matrC_CSR);
	if (status != 0) {
		cout << "matrix_mult_mkl failed" << endl;
		return -1;
	}
	MyTimer::SetEndTime();
	cout << "matrix multiplied using \"matrix_mult_mkl\" succesfully in " << MyTimer::GetDifferenceMs() << "ms" << endl;

	cout << "Matrix C:" << endl;
	cout << "N: " << matrC_CSR.N << ", M: " << matrC_CSR.M << ", nz: " << matrC_CSR.row_id[matrC_CSR.N] << endl;

	//mkl_sparse_destroy(matrA);
	//mkl_sparse_destroy(matrB);
	//mkl_sparse_destroy(matrC);
	//delete_CSR(&matrA_CSR);
	//delete_CSR(&matrB_CSR);
	//delete_CSR(&matrC_CSR);
	// how destoy correct ?

	transpose_this_CSR(&matrB_CSR); // !!! transpose matrix B !!!

	matrix_CSR matrC_CSR_naive;
	cout << "C_naive = A * B" << endl;
	MyTimer::SetStartTime();
	status = matrix_mult_naive(&matrA_CSR, &matrB_CSR, &matrC_CSR_naive);
	if (status != 0) {
		cout << "matrix_mult_naive failed" << endl;
		return -1;
	}
	MyTimer::SetEndTime();
	cout << "matrix multiplied using \"matrix_mult_naive\" succesfully in " << MyTimer::GetDifferenceMs() << "ms" << endl;

	cout << "Matrix C_naive:" << endl;
	cout << "N: " << matrC_CSR_naive.N << ", M: " << matrC_CSR_naive.M << ", nz: " << matrC_CSR_naive.row_id[matrC_CSR_naive.N] << endl;
	
	matrix_CSR matrC_CSR_copy;
	create_transposed_CSR(&matrC_CSR, &matrC_CSR_copy);
	transpose_this_CSR(&matrC_CSR_copy);
	transpose_this_CSR(&matrC_CSR_naive);
	transpose_this_CSR(&matrC_CSR_naive);
	if (matrix_compare(&matrC_CSR_copy, &matrC_CSR_naive)) {
		cout << "(!!!) matrC_CSR and matrC_CSR_naive are equal" << endl;
	}
	else {
		cout << "(!!!) Error: matrC_CSR != matrC_CSR_naive" << endl;
	}

	transpose_this_CSR(&matrA_CSR); // !!! matrA_CSR shoul be sorted !!! (matrB_CSR is already transposed and sorted)
	transpose_this_CSR(&matrA_CSR);

	matrix_CSR matrC_CSR_naive_2;
	cout << "C_naive_2 = A * B" << endl;
	MyTimer::SetStartTime();
	status = matrix_mult_naive_2(&matrA_CSR, &matrB_CSR, &matrC_CSR_naive_2);
	if (status != 0) {
		cout << "matrix_mult_naive failed" << endl;
		return -1;
	}
	MyTimer::SetEndTime();
	cout << "matrix multiplied using \"matrix_mult_naive_2\" succesfully in " << MyTimer::GetDifferenceMs() << "ms" << endl;

	cout << "Matrix C_naive:" << endl;
	cout << "N: " << matrC_CSR_naive_2.N << ", M: " << matrC_CSR_naive_2.M << ", nz: " << matrC_CSR_naive_2.row_id[matrC_CSR_naive_2.N] << endl;

	transpose_this_CSR(&matrC_CSR_naive_2);
	transpose_this_CSR(&matrC_CSR_naive_2);
	if (matrix_compare(&matrC_CSR_copy, &matrC_CSR_naive_2)) {
		cout << "(!!!) matrC_CSR and matrC_CSR_naive_2 are equal" << endl;
	}
	else {
		cout << "(!!!) Error: matrC_CSR != matrC_CSR_naive_2" << endl;
	}

	return 0;
}
