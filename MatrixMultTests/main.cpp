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

constexpr int TEST_NUM = 10;

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

int matrix_mult_TF(matrix_CSR* A_csr, matrix_CSR* B_csr, matrix_CSR* Res_csr, int delta) {
	int res = matrix_mult_TF_first_phase(A_csr, B_csr, Res_csr, delta);
	if (res != 0) return res;
	res = matrix_mult_TF_second_phase(A_csr, B_csr, Res_csr, delta);
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

	//my_test("naive_1_naiveomp", matrA_CSR, matrB_CSR, matrix_mult_naive_1_naiveomp, "mult_mkl", matrC_CSR);
	//my_test("naive_2_naiveomp", matrA_CSR, matrB_CSR, matrix_mult_naive_2_naiveomp, "mult_mkl", matrC_CSR);
	my_test("matrix_mult_naiveomp", matrA_CSR, matrB_CSR, matrix_mult_naiveomp, "mult_mkl", matrC_CSR);

	for (int delta : {50, 100, 250, 500}) {
		//my_test("naive_1_queueomp", matrA_CSR, matrB_CSR, matrix_mult_naive_1_queueomp, delta, "mult_mkl", matrC_CSR);
		//my_test("naive_2_queueomp", matrA_CSR, matrB_CSR, matrix_mult_naive_2_queueomp, delta, "mult_mkl", matrC_CSR);
		my_test("matrix_mult_queueomp", matrA_CSR, matrB_CSR, matrix_mult_queueomp, delta, "mult_mkl", matrC_CSR);

		my_test("matrix_mult_TF", matrA_CSR, matrB_CSR, matrix_mult_TF, delta, "mult_mkl", matrC_CSR);
	}

	delete_CSR(&matrC_CSR);

	cout << endl;
}

int main() {
	{
		cout << "-------------- small tests\n\n";
		matrix_COO a1_coo; read_matrix_MTX(R"(D:\source\coursework\sample_matrices\littlemult1.mtx)", &a1_coo);
		print_COO_as_dense(a1_coo, "a1_coo");
		cout << "a1_coo in COO:\n";
		cout << "I:\n";
		for (int i = 0; i < a1_coo.nz; i++) {
			cout << a1_coo.I[i] << " ";
		}
		cout << "\n";
		cout << "J:\n";
		for (int i = 0; i < a1_coo.nz; i++) {
			cout << a1_coo.J[i] << " ";
		}
		cout << "\n";
		cout << "val:\n";
		for (int i = 0; i < a1_coo.nz; i++) {
			cout << a1_coo.val[i] << " ";
		}
		cout << "\n";
		matrix_COO a2_coo; read_matrix_MTX(R"(D:\source\coursework\sample_matrices\littlemult2.mtx)", &a2_coo);
		print_COO_as_dense(a2_coo, "a2_coo");

		matrix_CSR a1; convert_COO_to_CSR(&a1_coo, &a1);
		matrix_CSR a2; convert_COO_to_CSR(&a2_coo, &a2);

		cout << "\na1 in CSR:\n";
		print_CSR(a1);

		matrix_CSR res; //matrix_mult_mkl(&a1, &a2, &res);
		transpose_this_CSR(&a2);
		matrix_mult_naive_1(&a1, &a2, &res);
		cout << "\n";
		print_CSR_as_dense(res, "res = a1 * a2");


		cout << "\n\n";
		matrix_CSR c1; generate_uniform_square_sparse_matrix(10, 3, c1);
		print_CSR_as_dense(c1, "c1, uniform");
		matrix_CSR c2; generate_nonuniform_square_sparse_matrix(10, 1, 4, c2);
		print_CSR_as_dense(c2, "c2, nonuniform");

		cout << "-------------- end of small tests\n\n";
	}

	{
		cout << "-------------------------------------------" << endl;
		cout << "test1 C = A * B" << endl;
		cout << "A: 10'000 x 10'000, 15 non-empty elements in every row" << endl;
		cout << "B: 10'000 x 10'000, 15 non-empty elements in every row" << endl;
		cout << "-------------------------------------------" << endl;
		cout << endl;

		matrix_CSR matrA_CSR;
		generate_uniform_square_sparse_matrix(10'000, 15, matrA_CSR);
		cout << "Matrix A generated successfully. N: " << matrA_CSR.N << ", M: " << matrA_CSR.M << ", nz: " << matrA_CSR.row_id[matrA_CSR.N] << endl;
		matrix_CSR matrB_CSR;
		generate_uniform_square_sparse_matrix(10'000, 15, matrB_CSR);
		cout << "Matrix B generated successfully. N: " << matrB_CSR.N << ", M: " << matrB_CSR.M << ", nz: " << matrB_CSR.row_id[matrB_CSR.N] << endl;
		cout << endl;

		basic_test(matrA_CSR, matrB_CSR);

		delete_CSR(&matrA_CSR);
		delete_CSR(&matrB_CSR);
	}

	{
		cout << "-------------------------------------------" << endl;
		cout << "test2 C = A * B" << endl;
		cout << "A: 10'000 x 10'000, from 30 to 1 non-elements in rows" << endl;
		cout << "B: 10'000 x 10'000, 15 non-empty elements in every row" << endl;
		cout << "-------------------------------------------" << endl;
		cout << endl;

		matrix_CSR matrA_CSR;
		generate_nonuniform_square_sparse_matrix(10'000, 30, 1, matrA_CSR);
		cout << "Matrix A generated successfully. N: " << matrA_CSR.N << ", M: " << matrA_CSR.M << ", nz: " << matrA_CSR.row_id[matrA_CSR.N] << endl;
		matrix_CSR matrB_CSR;
		generate_uniform_square_sparse_matrix(10'000, 15, matrB_CSR);
		cout << "Matrix B generated successfully" << endl;
		cout << "Matrix B generated successfully. N: " << matrB_CSR.N << ", M: " << matrB_CSR.M << ", nz: " << matrB_CSR.row_id[matrB_CSR.N] << endl;
		cout << endl;

		basic_test(matrA_CSR, matrB_CSR);

		delete_CSR(&matrA_CSR);
		delete_CSR(&matrB_CSR);
	}

	{
		cout << "-------------------------------------------" << endl;
		cout << "test3 C = A * B" << endl;
		cout << "A: 20'000 x 20'000, 25 non-empty elements in every row" << endl;
		cout << "B: 20'000 x 20'000, 25 non-empty elements in every row" << endl;
		cout << "-------------------------------------------" << endl;
		cout << endl;

		matrix_CSR matrA_CSR;
		generate_uniform_square_sparse_matrix(20'000, 25, matrA_CSR);
		cout << "Matrix A generated successfully. N: " << matrA_CSR.N << ", M: " << matrA_CSR.M << ", nz: " << matrA_CSR.row_id[matrA_CSR.N] << endl;
		matrix_CSR matrB_CSR;
		generate_uniform_square_sparse_matrix(20'000, 25, matrB_CSR);
		cout << "Matrix B generated successfully. N: " << matrB_CSR.N << ", M: " << matrB_CSR.M << ", nz: " << matrB_CSR.row_id[matrB_CSR.N] << endl;
		cout << endl;

		basic_test(matrA_CSR, matrB_CSR);

		delete_CSR(&matrA_CSR);
		delete_CSR(&matrB_CSR);
	}

	{
		cout << "-------------------------------------------" << endl;
		cout << "test4 C = A * B" << endl;
		cout << "A: 20'000 x 20'000, from 50 to 1 non-elements in rows" << endl;
		cout << "B: 20'000 x 20'000, 25 non-empty elements in every row" << endl;
		cout << "-------------------------------------------" << endl;
		cout << endl;

		matrix_CSR matrA_CSR;
		generate_nonuniform_square_sparse_matrix(20'000, 50, 1, matrA_CSR);
		cout << "Matrix A generated successfully. N: " << matrA_CSR.N << ", M: " << matrA_CSR.M << ", nz: " << matrA_CSR.row_id[matrA_CSR.N] << endl;
		matrix_CSR matrB_CSR;
		generate_uniform_square_sparse_matrix(20'000, 25, matrB_CSR);
		cout << "Matrix B generated successfully" << endl;
		cout << "Matrix B generated successfully. N: " << matrB_CSR.N << ", M: " << matrB_CSR.M << ", nz: " << matrB_CSR.row_id[matrB_CSR.N] << endl;
		cout << endl;

		basic_test(matrA_CSR, matrB_CSR);

		delete_CSR(&matrA_CSR);
		delete_CSR(&matrB_CSR);
	}

	return 0;
}
