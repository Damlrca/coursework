#ifndef MATRIX_IO_H
#define MATRIX_IO_H

int read_matrix_MTX(const char* fname, int* M_, int* N_, int* nz_, double** val_, int** I_, int** J_);

#endif // !MATRIX_IO_H
