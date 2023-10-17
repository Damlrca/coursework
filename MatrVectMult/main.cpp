#include <iostream>
extern "C" {
	#include "mmio.h"
}
#include <random>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <omp.h>
using namespace std;
using namespace std::chrono;
using myclock = chrono::system_clock;
myclock::time_point start_time, end_time;

const char* filename = "D:\\source\\coursework\\Freescale1.mtx";

void generate_CSR(int N, int M, int nz, double* val, int* I, int* J, int** _row_id, int** _col, double** _value) {
	int* row_id = new int[N + 1];
	int* col = new int[nz];
	double* value = new double[nz];

	int* cnt = new int[N] {};
	for (int i = 0; i < nz; i++) {
		cnt[I[i]]++;
	}
	int temp_id = 0;
	for (int i = 0; i < N; i++) {
		row_id[i] = temp_id;
		temp_id += cnt[i];
	}
	row_id[N] = temp_id;
	if (row_id[N] != nz) {
		cout << "error: row_id[N] != nz" << endl;
	}
	
	for (int i = 0; i < nz; i++) {
		col[row_id[I[i]]] = J[i];
		value[row_id[I[i]]] = val[i];
		++row_id[I[i]];
	}

	for (int i = 0; i < N; i++) {
		row_id[i] -= cnt[i];
	}

	delete[] cnt;

	*_row_id = row_id;
	*_col = col;
	*_value = value;
}

int main() {
	// read matrix
	start_time = myclock::now();
	int N, M, nz;
	int* I;
	int* J;
	double* val;
	if (mm_read_unsymmetric_sparse(filename, &N, &M, &nz, &val, &I, &J) != 0) {
		cout << "error in reading matrix" << endl;
		return -1;
	}
	end_time = myclock::now();
	cout << "matrix read succsessfully" << " in " << duration_cast<milliseconds>(end_time - start_time).count() << "ms" << endl;
	cout << "N: " << N << ", M: " << M << ", nz: " << nz << endl;

	// Generate vector V
	double* V = new double[M];

	std::random_device rd;
	std::mt19937_64 gen(rd());
	std::uniform_real_distribution<double> rand(-100.0, 100.0);
	for (int i = 0; i < M; i++) {
		V[i] = rand(gen);
	}

	// CSR format
	int* row_id;
	int* col;
	double* value;
	generate_CSR(N, M, nz, val, I, J, &row_id, &col, &value);

	// one thread, direct
	start_time = myclock::now();
	double* Result_direct = new double[N] {};
	for (int i = 0; i < nz; i++) {
		Result_direct[I[i]] += val[i] * V[J[i]];
	}
	end_time = myclock::now();
	cout << "one thread, direct calculated" << " in " << duration_cast<milliseconds>(end_time - start_time).count() << "ms" << endl;

	// one thread, CSR
	start_time = myclock::now();
	double* Result_CSR = new double[N] {};
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
	double Max_diff = abs(Result_direct[0] - Result_CSR[0]);
	for (int i = 0; i < N; i++)
		Max_diff = max(Max_diff, abs(Result_direct[i] - Result_CSR[i]));
	cout << "Max_diff (compare direct, CSR): " << Max_diff << endl;

	// some threads, easy openmp
	start_time = myclock::now();
	double* Result_easyOMP = new double[N] {};
	#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		int a = row_id[i];
		int b = row_id[i + 1];
		while (a != b) {
			Result_CSR[i] += value[a] * V[col[a]];
			++a;
		}
	}
	end_time = myclock::now();
	cout << "some threads, easy openmp calculated" << " in " << duration_cast<milliseconds>(end_time - start_time).count() << "ms" << endl;

	delete[] I;
	delete[] J;
	delete[] val;

	delete[] V;

	delete[] row_id;
	delete[] col;
	delete[] value;

	delete[] Result_direct;
	delete[] Result_CSR;
	delete[] Result_easyOMP;
	return 0;
}
