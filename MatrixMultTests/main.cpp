extern "C" {
#include "MatrixUtils.h"
#include "MatrixIO.h"
#include "MatrixMultMKL.h"
}
#include "MatrixGenerator.hpp"
#include "MatrixMult.h"
#include "MatrixMultTwoPhases.h"
#include "OtherUtils.hpp"
#include <iostream>
#include <limits>
#include <algorithm>

using std::cout;
using std::endl;

//const char* fileA = "D:\\source\\coursework\\sample_matrices\\multA.bin";
//const char* fileB = "D:\\source\\coursework\\sample_matrices\\multB.bin";

long long test_algo(const char* algo_name, int (*algo) (matrix_CSR*, matrix_CSR*, matrix_CSR*), matrix_CSR& A_csr, matrix_CSR& B_csr, matrix_CSR& C_csr) {
	MyTimer::SetStartTime();
	int status = algo(&A_csr, &B_csr, &C_csr);
	MyTimer::SetEndTime();
	if (status != 0) {
		cout << "(ERROR) : " << algo_name << " failed" << endl;
		exit(-1);
	}
	long long res = MyTimer::GetDifferenceMs();
	//cout << "(OK) : " << "matrix multiplied using \"" << algo_name << "\" succesfully in " << res << "ms" << endl;
	return res;
}

long long test_algo(const char* algo_name, int (*algo) (matrix_CSR*, matrix_CSR*, matrix_CSR*, int), int delta, matrix_CSR& A_csr, matrix_CSR& B_csr, matrix_CSR& C_csr) {
	MyTimer::SetStartTime();
	int status = algo(&A_csr, &B_csr, &C_csr, delta);
	MyTimer::SetEndTime();
	if (status != 0) {
		cout << "(ERROR) : " << algo_name << " failed" << endl;
		exit(-1);
	}
	long long res = MyTimer::GetDifferenceMs();
	//cout << "(OK) : " << "matrix multiplied using \"" << algo_name << "\" succesfully in " << res << "ms" << endl;
	return res;
}

void comp_algo(const char* algo_name_1, const char* algo_name_2, matrix_CSR& C_csr_1, matrix_CSR& C_csr_2) {
	if (matrix_compare(&C_csr_1, &C_csr_2)) {
		//cout << "(OK) : " << algo_name_1 << " == " << algo_name_2 << endl;
	}
	else {
		cout << "(ERROR) : " << algo_name_1 << " != " << algo_name_2 << endl;
		exit(-1);
	}
}

constexpr int TEST_NUM = 1;

void my_test(const char* algo_name, matrix_CSR& matrA_CSR, matrix_CSR& matrB_CSR, int (*algo) (matrix_CSR*, matrix_CSR*, matrix_CSR*),
	const char* algo_base_name, matrix_CSR& matrC_CSR_base) {
	matrix_CSR matrC_csr;
	long long best_ans = std::numeric_limits<long long>::max();
	for (int i = 0; i < TEST_NUM; i++) {
		long long t = test_algo(algo_name, algo, matrA_CSR, matrB_CSR, matrC_csr);
		best_ans = std::min(best_ans, t);
		transpose_this_CSR(&matrC_csr);
		transpose_this_CSR(&matrC_csr);
		comp_algo(algo_base_name, algo_name, matrC_CSR_base, matrC_csr);
		delete_CSR(&matrC_csr);
	}
	cout << "algo: " << algo_name << " : " << best_ans << "ms" << endl;
}

void my_test(const char* algo_name, matrix_CSR& matrA_CSR, matrix_CSR& matrB_CSR, int (*algo) (matrix_CSR*, matrix_CSR*, matrix_CSR*, int), int delta,
	const char* algo_base_name, matrix_CSR& matrC_CSR_base) {
	matrix_CSR matrC_csr;
	long long best_ans = std::numeric_limits<long long>::max();
	for (int i = 0; i < TEST_NUM; i++) {
		long long t = test_algo(algo_name, algo, delta, matrA_CSR, matrB_CSR, matrC_csr);
		best_ans = std::min(best_ans, t);
		transpose_this_CSR(&matrC_csr);
		transpose_this_CSR(&matrC_csr);
		comp_algo(algo_base_name, algo_name, matrC_CSR_base, matrC_csr);
		delete_CSR(&matrC_csr);
	}
	cout << "algo: " << algo_name << " delta: " << delta << " : " << best_ans << "ms" << endl;
}

int matrix_mult_TF_queueomp(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr, int delta) {
	int res = matrix_mult_TF_queueomp_first_phase(A_csr, B_csr, Res_csr, delta);
	if (res != 0) return res;
	res = matrix_mult_TF_queueomp_second_phase(A_csr, B_csr, Res_csr, delta);
	return res;
}

void basic_test(matrix_CSR& matrA_CSR, matrix_CSR& matrB_CSR) {
	matrix_CSR matrC_CSR;
	int basic_res = test_algo("matrix_mult_mkl", matrix_mult_mkl, matrA_CSR, matrB_CSR, matrC_CSR);
	transpose_this_CSR(&matrC_CSR);
	transpose_this_CSR(&matrC_CSR);
	cout << "Matrix C: N: " << matrC_CSR.N << ", M: " << matrC_CSR.M << ", nz: " << matrC_CSR.row_id[matrC_CSR.N] << endl;
	cout << endl;
	cout << "algo: mult_mkl : " << basic_res << "ms" << endl;

	transpose_this_CSR(&matrA_CSR); // !!! matrA_CSR should be sorted
	transpose_this_CSR(&matrA_CSR);
	transpose_this_CSR(&matrB_CSR); // !!! transpose matrix B !!!

	my_test("naive_1_naiveomp", matrA_CSR, matrB_CSR, matrix_mult_naive_1_naiveomp, "mult_mkl", matrC_CSR);
	my_test("naive_2_naiveomp", matrA_CSR, matrB_CSR, matrix_mult_naive_2_naiveomp, "mult_mkl", matrC_CSR);
	my_test("naive_3_naiveomp", matrA_CSR, matrB_CSR, matrix_mult_naive_3_naiveomp, "mult_mkl", matrC_CSR);

	for (int delta : {50, 100, 250, 500}) {
		my_test("naive_1_queueomp", matrA_CSR, matrB_CSR, matrix_mult_naive_1_queueomp, delta, "mult_mkl", matrC_CSR);
		my_test("naive_2_queueomp", matrA_CSR, matrB_CSR, matrix_mult_naive_2_queueomp, delta, "mult_mkl", matrC_CSR);
		my_test("naive_3_queueomp", matrA_CSR, matrB_CSR, matrix_mult_naive_3_queueomp, delta, "mult_mkl", matrC_CSR);

		my_test("TF_queueomp", matrA_CSR, matrB_CSR, matrix_mult_TF_queueomp, delta, "mult_mkl", matrC_CSR);
	}

	delete_CSR(&matrC_CSR);

	cout << endl;
}

int main() {
	{
		cout << "-------------------------------------------" << endl;
		cout << "test1 C = A * B" << endl;
		cout << "A: 10'000 x 10'000, 15 non-empty elements in every row" << endl;
		cout << "B: 10'000 x 10'000, 15 non-empty elements in every row" << endl;
		cout << "-------------------------------------------" << endl;
		cout << endl;

		matrix_CSR matrA_CSR;
		generate_uniform_square_sparse_matrix_CSR(10'000, 15, matrA_CSR);
		cout << "Matrix A generated successfully. N: " << matrA_CSR.N << ", M: " << matrA_CSR.M << ", nz: " << matrA_CSR.row_id[matrA_CSR.N] << endl;
		matrix_CSR matrB_CSR;
		generate_uniform_square_sparse_matrix_CSR(10'000, 15, matrB_CSR);
		cout << "Matrix B generated successfully. N: " << matrB_CSR.N << ", M: " << matrB_CSR.M << ", nz: " << matrB_CSR.row_id[matrB_CSR.N] << endl;
		cout << endl;

		basic_test(matrA_CSR, matrB_CSR);

		delete_CSR(&matrA_CSR);
		delete_CSR(&matrB_CSR);
	}

	{
		cout << "-------------------------------------------" << endl;
		cout << "test1 C = A * B" << endl;
		cout << "A: 10'000 x 10'000, from 30 to 1 non-elements in rows" << endl;
		cout << "B: 10'000 x 10'000, 15 non-empty elements in every row" << endl;
		cout << "-------------------------------------------" << endl;
		cout << endl;

		matrix_CSR matrA_CSR;
		generate_nonuniform_square_sparse_matrix_CSR(10'000, 30, 1, matrA_CSR);
		cout << "Matrix A generated successfully. N: " << matrA_CSR.N << ", M: " << matrA_CSR.M << ", nz: " << matrA_CSR.row_id[matrA_CSR.N] << endl;
		matrix_CSR matrB_CSR;
		generate_uniform_square_sparse_matrix_CSR(10'000, 15, matrB_CSR);
		cout << "Matrix B generated successfully" << endl;
		cout << "Matrix B generated successfully. N: " << matrB_CSR.N << ", M: " << matrB_CSR.M << ", nz: " << matrB_CSR.row_id[matrB_CSR.N] << endl;
		cout << endl;

		basic_test(matrA_CSR, matrB_CSR);

		delete_CSR(&matrA_CSR);
		delete_CSR(&matrB_CSR);
	}

	{
		cout << "-------------------------------------------" << endl;
		cout << "test1 C = A * B" << endl;
		cout << "A: 20'000 x 20'000, 25 non-empty elements in every row" << endl;
		cout << "B: 20'000 x 20'000, 25 non-empty elements in every row" << endl;
		cout << "-------------------------------------------" << endl;
		cout << endl;

		matrix_CSR matrA_CSR;
		generate_uniform_square_sparse_matrix_CSR(20'000, 25, matrA_CSR);
		cout << "Matrix A generated successfully. N: " << matrA_CSR.N << ", M: " << matrA_CSR.M << ", nz: " << matrA_CSR.row_id[matrA_CSR.N] << endl;
		matrix_CSR matrB_CSR;
		generate_uniform_square_sparse_matrix_CSR(20'000, 25, matrB_CSR);
		cout << "Matrix B generated successfully. N: " << matrB_CSR.N << ", M: " << matrB_CSR.M << ", nz: " << matrB_CSR.row_id[matrB_CSR.N] << endl;
		cout << endl;

		basic_test(matrA_CSR, matrB_CSR);

		delete_CSR(&matrA_CSR);
		delete_CSR(&matrB_CSR);
	}

	{
		cout << "-------------------------------------------" << endl;
		cout << "test1 C = A * B" << endl;
		cout << "A: 20'000 x 20'000, from 50 to 1 non-elements in rows" << endl;
		cout << "B: 20'000 x 20'000, 25 non-empty elements in every row" << endl;
		cout << "-------------------------------------------" << endl;
		cout << endl;

		matrix_CSR matrA_CSR;
		generate_nonuniform_square_sparse_matrix_CSR(20'000, 50, 1, matrA_CSR);
		cout << "Matrix A generated successfully. N: " << matrA_CSR.N << ", M: " << matrA_CSR.M << ", nz: " << matrA_CSR.row_id[matrA_CSR.N] << endl;
		matrix_CSR matrB_CSR;
		generate_uniform_square_sparse_matrix_CSR(20'000, 25, matrB_CSR);
		cout << "Matrix B generated successfully" << endl;
		cout << "Matrix B generated successfully. N: " << matrB_CSR.N << ", M: " << matrB_CSR.M << ", nz: " << matrB_CSR.row_id[matrB_CSR.N] << endl;
		cout << endl;

		basic_test(matrA_CSR, matrB_CSR);

		delete_CSR(&matrA_CSR);
		delete_CSR(&matrB_CSR);
	}

	return 0;
}
