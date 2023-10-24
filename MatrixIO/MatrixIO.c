#include "MatrixIO.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "mmio.h"

int read_matrix_MTX(const char* fname, int* M_, int* N_, int* nz_, double** val_, int** I_, int** J_) {
    // based on mm_read_unsymmetric_sparse

    FILE* f;
    MM_typecode matcode;
    int M, N, nz;
    int i;
    double* val;
    int* I;
    int* J;

    if ((f = fopen(fname, "r")) == NULL) {
        printf("read_matrix_MTX: Failed to read file %s\n", fname);
        return -1;
    }

    if (mm_read_banner(f, &matcode) != 0) {
        printf("read_matrix_MTX: Could not process Matrix Market banner in file [%s]\n", fname);
        return -1;
    }

    if (!(mm_is_real(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode))) {
        printf("read_matrix_MTX: This application does not support Matrix Market type: [%s]\n",
            mm_typecode_to_str(matcode));
        return -1;
    }

    if (mm_read_mtx_crd_size(f, &M, &N, &nz) != 0) {
        printf("read_matrix_MTX: Could not parse matrix size.\n");
        return -1;
    }

    if (mm_is_general(matcode)) {
        I = (int*)malloc(nz * sizeof(int));
        J = (int*)malloc(nz * sizeof(int));
        val = (double*)malloc(nz * sizeof(double));

        for (i = 0; i < nz; ++i) {
            if (fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]) != 3) {
                printf("read_matrix_MTX: Failed to read matrix data\n");
                return -1;
            }
            --I[i];
            --J[i];
        }

        *M_ = M;
        *N_ = N;
        *nz_ = nz;

        *val_ = val;
        *I_ = I;
        *J_ = J;
    }
    else { // matrix is symmetric
        int nnz = 0;
        int max_nz = nz * 2;

        I = (int*)malloc(max_nz * sizeof(int));
        J = (int*)malloc(max_nz * sizeof(int));
        val = (double*)malloc(max_nz * sizeof(double));

        for (i = 0; i < nz; ++i) {
            if (fscanf(f, "%d %d %lg\n", &I[nnz], &J[nnz], &val[nnz]) != 3) {
                printf("read_matrix_MTX: Failed to read matrix data\n");
                return -1;
            }
            --I[nnz];
            --J[nnz];
            ++nnz;
            if (I[nnz - 1] != J[nnz - 1]) {
                I[nnz] = J[nnz - 1];
                J[nnz] = I[nnz - 1];
                val[nnz] = val[nnz - 1];
                ++nnz;
            }
        }

        *M_ = M;
        *N_ = N;
        *nz_ = nnz;

        *val_ = val;
        *I_ = I;
        *J_ = J;
    }

    fclose(f);
    return 0;
}
