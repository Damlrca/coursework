#include <iostream>
extern "C" {
#include "MatrixIO.h"
#include "MatrixUtils.h"
}
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <random>
#include <algorithm>
#include <chrono>
//#include <omp.h>
using namespace std;
using namespace std::chrono;
using myclock = chrono::system_clock;
myclock::time_point start_time, end_time;

//const char* filename = "D:\\source\\coursework\\Freescale1.mtx";
const char* filename = "D:\\source\\coursework\\Freescale1.bin";

/*
void create_transposed(int N_A, int M_A, int* A_row_id, int* A_col, double* A_value,
	int* N_AT, int* M_AT, int** AT_row_id, int** AT_col, double** AT_value) {
	int nz = A_row_id[N_A];
	*N_AT = M_A;
	*M_AT = N_A;
	int* row_id = new int[*N_AT + 1];
	int* col = new int[nz];
	double* value = new double[nz];
	memset(row_id, 0, sizeof(int) * (*N_AT + 1));
	for (int i = 0; i < nz; i++) {
		row_id[col[i] + 1]++;
	}
	*AT_row_id = row_id;
	*AT_col = col;
	*AT_value = value;
}
*/

// TODO:
// сохранить матрицу в бинарный файл (для быстрого чтения, в частности)
// транспонирование матрицы
// сложение матриц
// умножение матриц

int main() {
	// read matrix
	start_time = myclock::now();
	int N, M, nz; // N - number of rows, M - number of columns
	int* I;
	int* J;
	double* val;
	//if (read_matrix_MTX(filename, &N, &M, &nz, &val, &I, &J) != 0) {
	if (read_matrix_BIN(filename, &N, &M, &nz, &val, &I, &J) != 0) {
		cout << "Failed to read matrix" << endl;
		return -1;
	}
	end_time = myclock::now();
	cout << "matrix read succsessfully" << " in " << duration_cast<milliseconds>(end_time - start_time).count() << "ms" << endl;
	cout << "N: " << N << ", M: " << M << ", nz: " << nz << endl;

	// CSR format
	int* row_id;
	int* col;
	double* value;
	COO_to_CSR(N, M, nz, val, I, J, &row_id, &col, &value);

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
	//for (int i = 0; i < 10; i++) {
	//	cout << Result_direct[i] << " " << Result_CSR[i] << endl;
	//}

	/*
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
	*/

	delete_matrix(&N, &M, &nz, &val, &I, &J);
	delete_CSR(&row_id, &col, &value);
	delete[] V;
	
	delete[] Result_direct;
	delete[] Result_CSR;
	//delete[] Result_easyOMP;
	return 0;
}
