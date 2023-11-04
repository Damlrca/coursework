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
	// small matrix 7x5 for tests
	const char* file_bin = "D:\\source\\coursework\\testmatrices\\littlematrix.bin";
	const char* file_mtx = "D:\\source\\coursework\\testmatrices\\littlematrix.mtx";

	// COO MTX, mU
	matrix_COO mtx1;
	read_matrix_MTX(file_mtx, &mtx1);
	cout << "Matrix " << file_mtx << endl;
	cout << "N: " << mtx1.N << ", M: " << mtx1.M << ", nz: " << mtx1.nz << endl;
	vector<vector<double>> mU(mtx1.N, vector<double>(mtx1.M));
	for (int i = 0; i < mtx1.nz; i++) {
		mU[mtx1.I[i]][mtx1.J[i]] = mtx1.val[i];
	}
	for (int i = 0; i < mtx1.N; i++) {
		for (int j = 0; j < mtx1.M; j++) {
			if (mU[i][j] >= 0) cout << " ";
			cout << mU[i][j] << " ";
		}
		cout << endl;
	}

	// COO BIN, mA
	matrix_COO mtx2;
	read_matrix_BIN(file_bin, &mtx2);
	cout << "Matrix " << file_bin << endl;
	cout << "N: " << mtx2.N << ", M: " << mtx2.M << ", nz: " << mtx2.nz << endl;
	vector<vector<double>> mA(mtx2.N, vector<double>(mtx2.M));
	for (int i = 0; i < mtx2.nz; i++) {
		mA[mtx2.I[i]][mtx2.J[i]] = mtx2.val[i];
	}
	for (int i = 0; i < mtx2.N; i++) {
		for (int j = 0; j < mtx2.M; j++) {
			if (mA[i][j] >= 0) cout << " ";
			cout << mA[i][j] << " ";
		}
		cout << endl;
	}

	// mU == mA ?
	for (int i = 0; i < mtx2.N; i++) {
		for (int j = 0; j < mtx2.M; j++) {
			if (mU[i][j] != mA[i][j]) {
				cout << "mU == mA: " << i << " " << j << " " << mU[i][j] << " " << mA[i][j] << endl;
			}
		}
	}

	// COO to CSR, mB
	matrix_CSR mtx3;
	COO_to_CSR(&mtx2, &mtx3);
	vector<vector<double>> mB(mtx3.N, vector<double>(mtx3.M));
	for (int i = 0; i < mtx3.N; i++) {
		int a = mtx3.row_id[i];
		int b = mtx3.row_id[i + 1];
		while (a != b) {
			mB[i][mtx3.col[a]] = mtx3.value[a];
			++a;
		}
	}
	cout << "Matrix COO_to_CSR:" << endl;
	cout << "N: " << mtx3.N << ", M: " << mtx3.M << endl;
	for (int i = 0; i < mtx3.N; i++) {
		for (int j = 0; j < mtx3.M; j++) {
			if (mB[i][j] >= 0) cout << " ";
			cout << mB[i][j] << " ";
		}
		cout << endl;
	}

	// mA == mB ?
	for (int i = 0; i < mtx3.N; i++) {
		for (int j = 0; j < mtx3.M; j++) {
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
	create_transposed(mtx2.N, mtx2.M, mtx3.row_id, mtx3.col, mtx3.value, &N_T, &M_T, &row_id_T, &col_T, &value_T);
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
	for (int i = 0; i < mtx2.N; i++) {
		for (int j = 0; j < mtx2.M; j++) {
			if (mA[i][j] != mT[j][i])
				cout << "mA != mT^T: " << i << " " << j << " " << mA[i][j] << " " << mT[j][i] << endl;
		}
	}

	delete_COO(&mtx1);
	delete_COO(&mtx2);
	delete_CSR(&mtx3);
	//_delete_CSR(&row_id_T, &col_T, &value_T);
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
	if (__read_matrix_BIN(filename, &N, &M, &nz, &val, &I, &J) != 0) {
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
	_COO_to_CSR(N, M, nz, val, I, J, &row_id, &col, &value);

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

	//_delete_COO(&val, &I, &J);
	//_delete_CSR(&row_id, &col, &value);
	N = 0; M = 0; nz = 0;
	delete[] V;

	delete[] Result_direct;
	delete[] Result_CSR;
	return 0;
}
