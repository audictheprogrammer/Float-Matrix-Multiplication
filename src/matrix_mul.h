#ifndef MATRIX_MUL_H
#define MATRIX_MUL_H

#include "matrix.h"

float modulo_naive(float a, int p);
void mp_naive(float** A, float** B, float** C, int n, int p);

void mp_ijk(float** A, float** B, float** C, int n);
void mp_kij(float** A, float** B, float** C, int n);
void mp_jki(float** A, float** B, float** C, int n);
void mp_ikj(float** A, float** B, float** C, int n);
void mp_jik(float** A, float** B, float** C, int n);
void mp_kji(float** A, float** B, float** C, int n);


#endif
