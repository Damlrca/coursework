#ifndef OTHER_UTILS_H
#define OTHER_UTILS_H

#include <iostream>
#include <vector>
#include <chrono>
extern "C" {
#include "MatrixUtils.h"
}

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
	static void PrintDifference() {
		std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << "ms";
	}
};

void print_dense(std::vector<std::vector<double>>& dense);

void print_COO_as_dense(const matrix_COO& mtx, const char* mtx_name);

void print_CSR(const matrix_CSR& mtx);

void print_CSR_as_dense(const matrix_CSR& mtx, const char* mtx_name);

void matrix_difference(matrix_CSR* mtx1, matrix_CSR* mtx2, matrix_CSR* res);

#endif // !OTHER_UTILS_H
