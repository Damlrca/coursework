extern "C" {
#include "MatrixUtils.h"
#include "MatrixIO.h"
#include "MatrixMultMKL.h"
}
#include "OtherUtils.hpp"
#include <iostream>

using std::cout;
using std::endl;

const char* fileA = "D:\\source\\coursework\\sample_matrices\\multA.bin";
const char* fileB = "D:\\source\\coursework\\sample_matrices\\multB.bin";

int main() {
	matrix_COO matrA_COO;
	read_matrix_BIN(fileA, &matrA_COO);
	matrix_COO matrB_COO;
	read_matrix_BIN(fileB, &matrB_COO);
	matrix_CSR matrA_CSR;
	convert_COO_to_CSR(&matrA_COO, &matrA_CSR);
	cout << "A (" << fileA << ")" << endl;
	cout << "Matrix A:" << endl;
	cout << "N: " << matrA_CSR.N << ", M: " << matrA_CSR.M << ", nz: " << matrA_CSR.row_id[matrA_CSR.N] << endl;
	matrix_CSR matrB_CSR;
	convert_COO_to_CSR(&matrB_COO, &matrB_CSR);
	cout << "B (" << fileA << ")" << endl;
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
	return 0;
}
