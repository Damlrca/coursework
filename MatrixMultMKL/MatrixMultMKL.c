#include "MatrixMultMKL.h"
#include "mkl_spblas.h"
#include <stdlib.h>
#include <string.h>

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
	matrix_CSR temp;
	mkl_sparse_d_export_csr(matrRes, &index_Res, &temp.N, &temp.M, &temp.row_id, &row_end_Res, &temp.col, &temp.value);

	Res_csr->N = temp.N;
	Res_csr->M = temp.M;
	Res_csr->row_id = (int*)malloc((Res_csr->N + 1) * sizeof(int));
	memcpy(Res_csr->row_id, temp.row_id, (Res_csr->N + 1) * sizeof(int));
	int nz = Res_csr->row_id[Res_csr->N];
	Res_csr->col = (int*)malloc(nz * sizeof(int));
	memcpy(Res_csr->col, temp.col, nz * sizeof(int));
	Res_csr->value = (double*)malloc(nz * sizeof(double));
	memcpy(Res_csr->value, temp.value, nz * sizeof(double));

	mkl_sparse_destroy(matrA);
	mkl_sparse_destroy(matrB);
	mkl_sparse_destroy(matrRes);
	return 0;
}
