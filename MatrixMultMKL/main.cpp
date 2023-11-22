extern "C" {
#include "MatrixUtils.h"
#include "MatrixIO.h"
}
#include "OtherUtils.hpp"
#include <iostream>
extern "C" {
#include "mkl_spblas.h"
}

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
	cout << matrA_CSR.N << " " << matrA_CSR.M << endl;
	cout << matrA_CSR.row_id[matrA_CSR.N] << endl;
	matrix_CSR matrB_CSR;
	convert_COO_to_CSR(&matrB_COO, &matrB_CSR);
	cout << matrB_CSR.N << " " << matrB_CSR.M << endl;
	cout << matrB_CSR.row_id[matrB_CSR.N] << endl;
	sparse_matrix_t matrA;
	mkl_sparse_d_create_csr(&matrA, SPARSE_INDEX_BASE_ZERO, matrA_CSR.N, matrA_CSR.M, matrA_CSR.row_id, matrA_CSR.row_id + 1, matrA_CSR.col, matrA_CSR.value);
	sparse_matrix_t matrB;
	mkl_sparse_d_create_csr(&matrB, SPARSE_INDEX_BASE_ZERO, matrB_CSR.N, matrB_CSR.M, matrB_CSR.row_id, matrB_CSR.row_id + 1, matrB_CSR.col, matrB_CSR.value);
	sparse_matrix_t matrC;
	sparse_status_t status;
	MyTimer::SetStartTime();
	status = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, matrA, matrB, &matrC);
	MyTimer::SetEndTime();
	if (status != SPARSE_STATUS_SUCCESS) {
		cout << "mkl_sparse_spmm failed" << endl;
		return -1;
	}
	cout << MyTimer::GetDifferenceMs() << "ms" << endl;
	matrix_CSR matrC_CSR;
	sparse_index_base_t index_C;
	int* row_end_C;
	mkl_sparse_d_export_csr(matrC, &index_C, &matrC_CSR.N, &matrC_CSR.M, &matrC_CSR.row_id, &row_end_C, &matrC_CSR.col, &matrC_CSR.value);
	cout << matrC_CSR.N << " " << matrC_CSR.M << endl;
	cout << matrC_CSR.row_id[matrC_CSR.N] << endl;

	mkl_sparse_destroy(matrA);
	mkl_sparse_destroy(matrB);
	mkl_sparse_destroy(matrC);
	//delete_CSR(&matrA_CSR);
	//delete_CSR(&matrB_CSR);
	//delete_CSR(&matrC_CSR);
	// how destoy correct ?
	return 0;
}
