#include <iostream>
extern "C" {
#include "MatrixIO.h"
}
#include <chrono>
using namespace std;
using namespace std::chrono;
using myclock = chrono::system_clock;
myclock::time_point start_time, end_time;

const char* filename_in = "D:\\source\\coursework\\sample_matrices\\littlematrix.mtx";
const char* filename_out = "D:\\source\\coursework\\sample_matrices\\littlematrix.bin";

int main() {
	// read MTX matrix
	start_time = myclock::now();
	matrix_COO mtx;
	if (read_matrix_MTX(filename_in, &mtx) != 0) {
		cout << "Failed to read matrix" << endl;
		return -1;
	}
	end_time = myclock::now();
	cout << "matrix read succsessfully" << " in " << duration_cast<milliseconds>(end_time - start_time).count() << "ms" << endl;
	cout << "N: " << mtx.N << ", M: " << mtx.M << ", nz: " << mtx.nz << endl;

	// write BIN watrix
	start_time = myclock::now();
	if (save_matrix_BIN(filename_out, &mtx) != 0) {
		cout << "Failed to save matrix" << endl;
		return -1;
	}
	end_time = myclock::now();
	cout << "matrix saved succsessfully" << " in " << duration_cast<milliseconds>(end_time - start_time).count() << "ms" << endl;
	return 0;
}
