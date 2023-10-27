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
#include <vector>
using namespace std;
using namespace std::chrono;
using myclock = chrono::system_clock;
myclock::time_point start_time, end_time;

//const char* filename = "D:\\source\\coursework\\Freescale1.mtx";
const char* filename = "D:\\source\\coursework\\Freescale1.bin";

void small_test() {
	cout << fixed; cout.precision(2);
	// small matrix 7x7 for tests
	const char* file_bin = "D:\\source\\coursework\\testmatrices\\littlematrix.bin";
	const char* file_mtx = "D:\\source\\coursework\\testmatrices\\littlematrix.mtx";

	// COO MTX, mU
	int N_mtx, M_mtx, nz_mtx;
	double* val_mtx;
	int* I_mtx;
	int* J_mtx;
	read_matrix_MTX(file_mtx, &N_mtx, &M_mtx, &nz_mtx, &val_mtx, &I_mtx, &J_mtx);
	cout << "Matrix" << " " << file_mtx << endl;
	cout << "N: " << N_mtx << ", M: " << M_mtx << ", nz: " << nz_mtx << endl;
	vector<vector<double>> mU(N_mtx, vector<double>(M_mtx));
	for (int i = 0; i < nz_mtx; i++) {
		mU[I_mtx[i]][J_mtx[i]] = val_mtx[i];
	}
	for (int i = 0; i < N_mtx; i++) {
		for (int j = 0; j < M_mtx; j++) {
			if (mU[i][j] >= 0) cout << " ";
			cout << mU[i][j] << " ";
		}
		cout << endl;
	}

	// COO BIN, mA
	int N, M, nz;
	double* val;
	int* I;
	int* J;
	read_matrix_BIN(file_bin, &N, &M, &nz, &val, &I, &J);
	cout << "Matrix" << " " << file_bin << endl;
	cout << "N: " << N << ", M: " << M << ", nz: " << nz << endl;
	vector<vector<double>> mA(N, vector<double>(M));
	for (int i = 0; i < nz; i++) {
		mA[I[i]][J[i]] = val[i];
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			if (mA[i][j] >= 0) cout << " ";
			cout << mA[i][j] << " ";
		}
		cout << endl;
	}

	// mU == mA ?
	for (int i = 0; i < N_mtx; i++) {
		for (int j = 0; j < M_mtx; j++) {
			if (mU[i][j] != mA[i][j]) {
				cout << "mU == mA: " << i << " " << j << " " << mU[i][j] << " " << mA[i][j] << endl;
			}
		}
	}

	// COO to CSR, mB
	int* row_id;
	int* col;
	double* value;
	COO_to_CSR(N, M, nz, val, I, J, &row_id, &col, &value);
	vector<vector<double>> mB(N, vector<double>(M));
	for (int i = 0; i < N; i++) {
		int a = row_id[i];
		int b = row_id[i + 1];
		while (a != b) {
			mB[i][col[a]] = value[a];
			++a;
		}
	}
	cout << "CSR:" << endl;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			if (mA[i][j] >= 0) cout << " ";
			cout << mA[i][j] << " ";
		}
		cout << endl;
	}

	// mA == mB ?
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			if (mA[i][j] != mB[i][j]) {
				cout << "mA != mB: " << i << " " << j << " " << mA[i][j] << " " << mB[i][j] << endl;
			}
		}
	}

	// CSR transposed
	int N_T, M_T;
	int* row_id_T;
	int* col_T;
	double* value_T;
	create_transposed(N, M, row_id, col, value, &N_T, &M_T, &row_id_T, &col_T, &value_T);
	cout << "transposed: ";
	cout << "N_T: " << N_T << ", M_T: " << M_T << ", nz_T: " << row_id_T[N_T] << endl;
	vector<vector<double>> mT(N_T, vector<double>(M_T));
	for (int i = 0; i < N_T; i++) {
		int a = row_id_T[i];
		int b = row_id_T[i + 1];
		while (a != b) {
			mT[i][col_T[a]] = value_T[a];
			++a;
		}
	}
	cout << "CSR transposed:" << endl;
	for (int i = 0; i < N_T; i++) {
		for (int j = 0; j < M_T; j++) {
			if (mT[i][j] >= 0) cout << " ";
			cout << mT[i][j] << " ";
		}
		cout << endl;
	}

	// mA == mT^T ?
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			if (mA[i][j] != mT[j][i])
				cout << "mA != mT^T: " << i << " " << j << " " << mA[i][j] << " " << mT[j][i] << endl;
		}
	}

	delete_COO(&val_mtx, &I_mtx, &J_mtx);
	delete_COO(&val, &I, &J);
	delete_CSR(&row_id, &col, &value);
	delete_CSR(&row_id_T, &col_T, &value_T);
};

int main() {
	small_test();
	return 0;

	// read matrix
	start_time = myclock::now();
	int N, M, nz; // N - number of rows, M - number of columns
	int* I;
	int* J;
	double* val;
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

	delete_COO(&val, &I, &J);
	delete_CSR(&row_id, &col, &value);
	N = 0; M = 0; nz = 0;
	delete[] V;

	delete[] Result_direct;
	delete[] Result_CSR;
	return 0;
}
