#include <iostream>
extern "C" {
#include "MatrixIO.h"
#include "MatrixUtils.h"
}
#include "OtherUtils.hpp"
#include <iostream>
#include <chrono>
using namespace std;

const char* filename_in = "D:\\source\\coursework\\sample_matrices\\littlematrix.mtx";
const char* filename_out = "D:\\source\\coursework\\sample_matrices\\littlematrix.bin";

int main() {
	// read MTX matrix
	MyTimer::SetStartTime();
	matrix_COO mtx;
	if (read_matrix_MTX(filename_in, &mtx) != 0) {
		cout << "Failed to read matrix" << endl;
		return -1;
	}
	MyTimer::SetEndTime();
	cout << "matrix read succsessfully in " << MyTimer::GetDifferenceMs() << "ms" << endl;
	cout << "N: " << mtx.N << ", M: " << mtx.M << ", nz: " << mtx.nz << endl;

	// write BIN watrix
	MyTimer::SetStartTime();
	if (save_matrix_BIN(filename_out, &mtx) != 0) {
		cout << "Failed to save matrix" << endl;
		return -1;
	}
	MyTimer::SetEndTime();
	cout << "matrix saved succsessfully in " << MyTimer::GetDifferenceMs() << "ms" << endl;
	return 0;
}
