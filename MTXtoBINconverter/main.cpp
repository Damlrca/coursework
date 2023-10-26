#include <iostream>
extern "C" {
#include "MatrixIO.h"
}
#include <chrono>
using namespace std;
using namespace std::chrono;
using myclock = chrono::system_clock;
myclock::time_point start_time, end_time;

const char* filename_in = "D:\\source\\coursework\\Freescale1.mtx";
const char* filename_out = "D:\\source\\coursework\\Freescale1.bin";

int main() {
	// read MTX matrix
	start_time = myclock::now();
	int N, M, nz; // N - number of rows, M - number of columns
	int* I;
	int* J;
	double* val;
	if (read_matrix_MTX(filename_in, &N, &M, &nz, &val, &I, &J) != 0) {
		cout << "Failed to read matrix" << endl;
		return -1;
	}
	end_time = myclock::now();
	cout << "matrix read succsessfully" << " in " << duration_cast<milliseconds>(end_time - start_time).count() << "ms" << endl;
	cout << "N: " << N << ", M: " << M << ", nz: " << nz << endl;

	// write BIN watrix
	start_time = myclock::now();
	if (save_matrix_BIN(filename_out, N, M, nz, val, I, J) != 0) {
		cout << "Failed to save matrix" << endl;
		return -1;
	}
	end_time = myclock::now();
	cout << "matrix saved succsessfully" << " in " << duration_cast<milliseconds>(end_time - start_time).count() << "ms" << endl;
	return 0;
}
