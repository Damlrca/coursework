extern "C" {
#include "MatrixIO.h"
#include "MatrixUtils.h"
}
#include "OtherUtils.hpp"
#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <random>
#include <algorithm>
#include <vector>
using namespace std;

void few_tests() {
	// small matrix 7x5 for tests
	const char* file_bin = "D:\\source\\coursework\\sample_matrices\\littlematrix.bin";
	const char* file_mtx = "D:\\source\\coursework\\sample_matrices\\littlematrix.mtx";

	// COO MTX (mtx1)
	matrix_COO mtx1;
	read_matrix_MTX(file_mtx, &mtx1);
	print_COO_as_dense(mtx1, "read COO mtx : mtx1");
	cout << "----------------------" << endl;

	// COO BIN (mtx2)
	matrix_COO mtx2;
	read_matrix_BIN(file_bin, &mtx2);
	print_COO_as_dense(mtx2, "read COO bin : mtx2");
	cout << "----------------------" << endl;

	// COO to CSR (mtx3)
	matrix_CSR mtx3;
	convert_COO_to_CSR(&mtx2, &mtx3);
	print_CSR_as_dense(mtx3, "convert COO to CSR : mtx3");
	cout << "----------------------" << endl;

	// CSR transposed (mtx4)
	matrix_CSR mtx4;
	create_transposed_CSR(&mtx3, &mtx4);
	print_CSR_as_dense(mtx4, "create transposed CSR : mtx4");
	cout << "----------------------" << endl;

	// transpose this (mtx4)
	transpose_this_CSR(&mtx4);
	print_CSR_as_dense(mtx4, "transpose this CSR : mtx4");
	cout << "----------------------" << endl;

	if (matrix_compare(&mtx3, &mtx4) == false) {
		cout << "something failed: mtx3 != mtx4  :(" << endl;
		exit(-1);
	}

	delete_COO(&mtx1);
	delete_COO(&mtx2);
	delete_CSR(&mtx3);
	delete_CSR(&mtx4);
};

void test_addition() {
	const char* file1 = "D:\\source\\coursework\\sample_matrices\\add1.mtx";
	const char* file2 = "D:\\source\\coursework\\sample_matrices\\add2.mtx";
	
	matrix_COO add1_COO;
	read_matrix_MTX(file1, &add1_COO);
	print_COO_as_dense(add1_COO, "read add1 matrix");
	cout << "----------------------" << endl;

	matrix_COO add2_COO;
	read_matrix_MTX(file2, &add2_COO);
	print_COO_as_dense(add2_COO, "read add2 matrix");
	cout << "----------------------" << endl;

	matrix_CSR add1;
	convert_COO_to_CSR(&add1_COO, &add1);
	
	cout << "add1 in CSR format" << endl;
	print_CSR(add1);
	cout << "----------------------" << endl;
	transpose_this_CSR(&add1);
	transpose_this_CSR(&add1);
	cout << "add1 in ORDERED CSR format" << endl;
	print_CSR(add1);
	cout << "----------------------" << endl;

	matrix_CSR add2;
	convert_COO_to_CSR(&add2_COO, &add2);
	cout << "add2 in CSR format" << endl;
	print_CSR(add2);
	cout << "----------------------" << endl;
	transpose_this_CSR(&add2);
	transpose_this_CSR(&add2);
	cout << "add2 in ORDERED CSR format" << endl;
	print_CSR(add2);
	cout << "----------------------" << endl;

	matrix_CSR res;
	matrix_addition(&add1, &add2, &res);
	print_CSR_as_dense(res, "result");
	cout << "----------------------" << endl;
	cout << "result in CSR format" << endl;
	print_CSR(res);
	cout << "----------------------" << endl;

	delete_COO(&add1_COO);
	delete_COO(&add2_COO);
	delete_CSR(&add1);
	delete_CSR(&add2);
	delete_CSR(&res);
}

void few_more_tests() {
	const char* filename_mtx = "D:\\source\\coursework\\sample_matrices\\Freescale1.mtx";
	const char* filename_bin = "D:\\source\\coursework\\sample_matrices\\Freescale1.bin";

	// read mtx matrix
	matrix_COO mtx1;
	cout << "(mtx1) read " << filename_mtx << endl;
	MyTimer::SetStartTime();
	if (read_matrix_MTX(filename_mtx, &mtx1) != 0) {
		cout << "failed to read " << filename_mtx << endl;
		exit(-1);
	};
	MyTimer::SetEndTime();
	cout << "N :" << mtx1.N << ", M: " << mtx1.M << ", nz: " << mtx1.nz << endl;
	cout << "matrix read successfully in " << MyTimer::GetDifferenceMs() << "ms" << endl;
	cout << "----------------------" << endl;
	
	// read bin matrix
	matrix_COO mtx2;
	cout << "(mtx2) read " << filename_bin << endl;
	MyTimer::SetStartTime();
	if (read_matrix_BIN(filename_bin, &mtx2) != 0) {
		cout << "failed to read" << " " << filename_bin << endl;
		exit(-1);
	};
	MyTimer::SetEndTime();
	cout << "N :" << mtx2.N << ", M: " << mtx2.M << ", nz: " << mtx2.nz << endl;
	cout << "matrix read successfully in " << MyTimer::GetDifferenceMs() << "ms" << endl;
	cout << "----------------------" << endl;

	// mtx COO to CSR
	cout << "(mtx1->mtx3) convert COO to CSR" << endl;
	matrix_CSR mtx3;
	MyTimer::SetStartTime();
	convert_COO_to_CSR(&mtx1, &mtx3);
	MyTimer::SetEndTime();
	cout << "N :" << mtx3.N << ", M: " << mtx3.M << ", nz: " << mtx3.row_id[mtx3.N] << endl;
	cout << "matrix converted successfully in " << MyTimer::GetDifferenceMs() << "ms" << endl;
	cout << "----------------------" << endl;

	// bin COO to CSR
	cout << "(mtx2->mtx4) convert COO to CSR" << endl;
	matrix_CSR mtx4;
	MyTimer::SetStartTime();
	convert_COO_to_CSR(&mtx2, &mtx4);
	MyTimer::SetEndTime();
	cout << "N :" << mtx4.N << ", M: " << mtx4.M << ", nz: " << mtx4.row_id[mtx4.N] << endl;
	cout << "matrix converted successfully in " << MyTimer::GetDifferenceMs() << "ms" << endl;
	cout << "----------------------" << endl;

	// compare matrices mtx3 and mtx4
	/* // to compare two matrices using matrix_compare function they should be in ORDERED CSR FORMAT
	cout << "comparing mtx3 and mtx4" << endl;
	bool result;
	MyTimer::SetStartTime();
	result = matrix_compare(&mtx3, &mtx4);
	MyTimer::SetEndTime();
	cout << "matrices compared in " << MyTimer::GetDifferenceMs() << "ms" << endl;
	if (result) {
		cout << "OK, mtx3 == mtx4" << endl;
	}
	else {
		cout << "something failed: mtx3 != mtx4  :(" << endl;
		exit(-1);
	}
	cout << "----------------------" << endl;
	*/

	// double transpose mtx4
	cout << "(mtx4) double transpose" << endl;
	MyTimer::SetStartTime();
	transpose_this_CSR(&mtx4);
	transpose_this_CSR(&mtx4);
	MyTimer::SetEndTime();
	cout << "N :" << mtx4.N << ", M: " << mtx4.M << ", nz: " << mtx4.row_id[mtx4.N] << endl;
	cout << "matrix double transposed successfully in " << MyTimer::GetDifferenceMs() << "ms" << endl;
	cout << "----------------------" << endl;

	// compare matrices mtx3 ans mtx4 after double transpose
	/* // to compare two matrices using matrix_compare function they should be in ORDERED CSR FORMAT
	cout << "comparing mtx3 and mtx4" << endl;
	MyTimer::SetStartTime();
	result = matrix_compare(&mtx3, &mtx4);
	MyTimer::SetEndTime();
	cout << "matrices compared in " << MyTimer::GetDifferenceMs() << "ms" << endl;
	if (result) {
		cout << "OK, mtx3 == mtx4" << endl;
	}
	else {
		cout << "something failed: mtx3 != mtx4  :(" << endl;
		exit(-1);
	}
	cout << "----------------------" << endl;
	*/

	delete_COO(&mtx1);
	delete_COO(&mtx2);
	delete_CSR(&mtx3);
	delete_CSR(&mtx4);
}

int main() {
	cout << fixed; cout.precision(2);
	cout << "Testing main functions" << endl;
	cout << "----------------------" << endl;
	few_tests();
	cout << "Testing addition" << endl;
	cout << "----------------------" << endl;
	test_addition();
	cout << "Testing some functions with larger matrices" << endl;
	cout << "----------------------" << endl;
	few_more_tests();
	return 0;
}
