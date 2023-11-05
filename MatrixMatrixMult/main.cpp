extern "C" {
#include "MatrixIO.h"
#include "MatrixUtils.h"
}
#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <random>
#include <algorithm>
#include <chrono>
#include <vector>
#include <iomanip>
using namespace std;

const char* filename_mtx = "D:\\source\\coursework\\Freescale1.mtx";
const char* filename_bin = "D:\\source\\coursework\\Freescale1.bin";

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
		cout << chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count() << "ms";
	}
};
MyTimer::myclock::time_point MyTimer::start_time;
MyTimer::myclock::time_point MyTimer::end_time;

void print_dense(vector<vector<double>>& dense) {
	for (int i = 0; i < dense.size(); i++) {
		for (int j = 0; j < dense[i].size(); j++)
			cout << setw(5) << dense[i][j] << " ";
		cout << endl;
	}
}

void print_COO_as_dense(const matrix_COO& mtx, const char* mtx_name) {
	cout << mtx_name << endl;
	cout << "N: " << mtx.N << ", M: " << mtx.M << ", nz: " << mtx.nz << endl;
	vector<vector<double>> dense(mtx.N, vector<double>(mtx.M));
	for (int i = 0; i < mtx.nz; i++) {
		dense[mtx.I[i]][mtx.J[i]] = mtx.val[i];
	}
	print_dense(dense);
}

void print_CSR_as_dense(const matrix_CSR& mtx, const char* mtx_name) {
	cout << mtx_name << endl;
	cout << "N: " << mtx.N << ", M: " << mtx.M << ", nz: " << mtx.row_id[mtx.N] << endl;
	vector<vector<double>> dense(mtx.N, vector<double>(mtx.M));
	for (int i = 0; i < mtx.N; i++) {
		int a = mtx.row_id[i];
		int b = mtx.row_id[i + 1];
		while (a != b) {
			dense[i][mtx.col[a]] = mtx.value[a];
			++a;
		}
	}
	print_dense(dense);
}

void few_tests() {
	cout << fixed; cout.precision(2);
	// small matrix 7x5 for tests
	const char* file_bin = "D:\\source\\coursework\\testmatrices\\littlematrix.bin";
	const char* file_mtx = "D:\\source\\coursework\\testmatrices\\littlematrix.mtx";

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
	COO_to_CSR(&mtx2, &mtx3);
	print_CSR_as_dense(mtx3, "convert COO to CSR : mtx3");
	cout << "----------------------" << endl;

	// CSR transposed (mtx4)
	matrix_CSR mtx4;
	create_transposed(&mtx3, &mtx4);
	print_CSR_as_dense(mtx4, "create transposed CSR : mtx4");
	cout << "----------------------" << endl;

	// transpose this (mtx4)
	transpose_this(&mtx4);
	print_CSR_as_dense(mtx4, "transpose this CSR : mtx4");
	cout << "----------------------" << endl;

	delete_COO(&mtx1);
	delete_COO(&mtx2);
	delete_CSR(&mtx3);
	delete_CSR(&mtx4);
};

void few_more_tests() {
	// read mtx matrix
	cout << "read " << filename_mtx << endl;
	MyTimer::SetStartTime();
	matrix_COO mtx1;
	if (read_matrix_MTX(filename_mtx, &mtx1) != 0) {
		cout << "failed to read" << " " << filename_mtx << endl;
		exit(-1);
	};
	MyTimer::SetEndTime();
	cout << "N :" << mtx1.N << ", M: " << mtx1.M << ", nz: " << mtx1.nz << endl;
	cout << "matrix read succsessfully in ";
	MyTimer::PrintDifference();
	cout << endl;
	cout << "----------------------" << endl;

	// read bin matrix
	cout << "read " << filename_bin << endl;
	MyTimer::SetStartTime();
	matrix_COO mtx2;
	if (read_matrix_BIN(filename_bin, &mtx2) != 0) {
		cout << "failed to read" << " " << filename_bin << endl;
		exit(-1);
	};
	MyTimer::SetEndTime();
	cout << "N :" << mtx2.N << ", M: " << mtx2.M << ", nz: " << mtx2.nz << endl;
	cout << "matrix read succsessfully in ";
	MyTimer::PrintDifference();
	cout << endl;
	cout << "----------------------" << endl;

	// COO to CSR
	cout << "convert COO to CSR" << endl;
	matrix_CSR mtx3;
	MyTimer::SetStartTime();
	COO_to_CSR(&mtx2, &mtx3);
	MyTimer::SetEndTime();
	cout << "N :" << mtx3.N << ", M: " << mtx3.M << ", nz: " << mtx3.row_id[mtx3.N] << endl;
	cout << "matrix converted succsessfully in ";
	MyTimer::PrintDifference();
	cout << endl;
	cout << "----------------------" << endl;

	delete_COO(&mtx1);
	delete_COO(&mtx2);
	delete_CSR(&mtx3);
}

int main() {
	cout << "Testing main functions" << endl;
	cout << "----------------------" << endl;
	few_tests();
	cout << "Testing some functions with larger matrices" << endl;
	cout << "----------------------" << endl;
	few_more_tests();
	return 0;
	/*
	// check matrix-vector multiplication

	// Generate vector V
	double* V = new double[M];
	std::random_device rd;
	std::mt19937_64 gen(rd());
	std::uniform_real_distribution<double> rand(-100.0, 100.0);
	for (int i = 0; i < M; i++)
		V[i] = rand(gen);

	// one thread, direct
	start_time = myclock::now();
	double* Result_direct = new double[N];
	memset(Result_direct, 0, sizeof(double) * N);
	for (int i = 0; i < nz; i++) {
		Result_direct[I[i]] += val[i] * V[J[i]];
	}
	end_time = myclock::now();
	cout << "one thread, direct calculated" << " in " << duration_cast<milliseconds>(end_time - start_time).count() << "ms" << endl;

	// one thread, CSR
	start_time = myclock::now();
	double* Result_CSR = new double[N];
	memset(Result_CSR, 0, sizeof(double) * N);
	for (int i = 0; i < N; i++) {
		int a = row_id[i];
		int b = row_id[i + 1];
		while (a != b) {
			Result_CSR[i] += value[a] * V[col[a]];
			++a;
		}
	}
	end_time = myclock::now();
	cout << "one thread, CSR calculated" << " in " << duration_cast<milliseconds>(end_time - start_time).count() << "ms" << endl;

	// compare direct and CSR
	double Max_diff = 0;
	for (int i = 0; i < N; i++)
		Max_diff = max(Max_diff, abs(Result_direct[i] - Result_CSR[i]));
	cout << "Max_diff in result vectors (direct, CSR): " << Max_diff << endl;

	N = 0; M = 0; nz = 0;
	delete[] V;

	delete[] Result_direct;
	delete[] Result_CSR;
	*/
	return 0;
}
