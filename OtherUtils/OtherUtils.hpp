#ifndef OTHER_UTILS_H
#define OTHER_UTILS_H

#include <iostream>
#include <vector>
#include <chrono>
extern "C" {
#include "MatrixUtils.h"
}

// class for measuring time
// MyTimer::SetStartTime();
// ... code ...
// MyTimer::SetEndTime();
// cout << MyTimer::GetDifferenceMs() << "ms" << endl;
class MyTimer {
	using myclock = std::chrono::system_clock;
	static myclock::time_point start_time;
	static myclock::time_point end_time;
public:
	static void SetStartTime() {
		start_time = myclock::now();
	}
	static void SetEndTime() {
		end_time = myclock::now();
	}
	static long long GetDifferenceMs() {
		return std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
	}
};

// print dense matrix
void print_dense(std::vector<std::vector<double>>& dense);

// print COO-matrix as dense matrix
void print_COO_as_dense(const matrix_COO& mtx, const char* mtx_name);

// print CSR-matrix in CSR format
void print_CSR(const matrix_CSR& mtx);

// print CSR-matrix as dense matrix
void print_CSR_as_dense(const matrix_CSR& mtx, const char* mtx_name);

// difference of two CSR-matrices
// matrices "mtx1" and "mtx2" should be ordered!
// matrix "res" should be empty!
void matrix_difference(matrix_CSR* mtx1, matrix_CSR* mtx2, matrix_CSR* res);

// returns true if matrices "mtx1" and "mtx2" are equal
// matrices "mtx1" and "mtx2" should be ordered!
bool matrix_compare(matrix_CSR* mtx1, matrix_CSR* mtx2);

#endif // !OTHER_UTILS_H
