#include "OtherUtils.hpp"
#include <iomanip>
using namespace std;

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

void print_CSR(const matrix_CSR& mtx) {
	int N = mtx.N;
	int M = mtx.M;
	int nz = mtx.row_id[N];
	cout << N << " " << M << " " << nz << endl;
	for (int i = 0; i < N + 1; i++) {
		cout << mtx.row_id[i] << " ";
	}
	cout << endl;
	for (int i = 0; i < nz; i++) {
		cout << mtx.col[i] << " ";
	}
	cout << endl;
	for (int i = 0; i < nz; i++) {
		cout << setw(5) << mtx.value[i] << " ";
	}
	cout << endl;
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

// difference of two matricies in ORDERED CSR format
void matrix_difference(matrix_CSR* mtx1, matrix_CSR* mtx2, matrix_CSR* res) {
	int nz = mtx2->row_id[mtx2->N];
	for (int i = 0; i < nz; i++)
		mtx2->value[i] = -mtx2->value[i];
	matrix_addition(mtx1, mtx2, res);
	for (int i = 0; i < nz; i++)
		mtx2->value[i] = -mtx2->value[i];
}
