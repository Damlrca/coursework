#include "MatrixMultMKL.h"
#include "mkl_spblas.h"

int matrix_mult_mkl(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr) {
	sparse_matrix_t matrA;
	mkl_sparse_d_create_csr(&matrA, SPARSE_INDEX_BASE_ZERO, A_csr->N, A_csr->M, A_csr->row_id, A_csr->row_id + 1, A_csr->col, A_csr->value);
	sparse_matrix_t matrB;
	mkl_sparse_d_create_csr(&matrB, SPARSE_INDEX_BASE_ZERO, B_csr->N, B_csr->M, B_csr->row_id, B_csr->row_id + 1, B_csr->col, B_csr->value);

	sparse_matrix_t matrRes;
	sparse_status_t status;
	status = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, matrA, matrB, &matrRes);
	if (status != SPARSE_STATUS_SUCCESS) {
		return -1;
	}

	sparse_index_base_t index_Res;
	int* row_end_Res;
	mkl_sparse_d_export_csr(matrRes, &index_Res, &Res_csr->N, &Res_csr->M, &Res_csr->row_id, &row_end_Res, &Res_csr->col, &Res_csr->value);

	mkl_sparse_destroy(matrA);
	mkl_sparse_destroy(matrB);
	//mkl_sparse_destroy(matrC);
	//delete_CSR(&matrA_CSR);
	//delete_CSR(&matrB_CSR);
	//delete_CSR(&matrC_CSR);
	// TODO : how destoy correct ?
	return 0;
}
