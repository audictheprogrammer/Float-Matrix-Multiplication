#ifndef MATRIX_MUL_H
#define MATRIX_MUL_H

#include "matrix.h"

double modulo_naive(double a, double p);
double modulo_SIMD1(double a, double p, double u);

void mp_naive(double** A, double** B, double** C, int n, double p);
void mp_SIMD1(double** A, double** B, double** C, int n, double p, double u);

void mp_ijk(double** A, double** B, double** C, int n);
void mp_kij(double** A, double** B, double** C, int n);
void mp_jki(double** A, double** B, double** C, int n);
void mp_ikj(double** A, double** B, double** C, int n);
void mp_jik(double** A, double** B, double** C, int n);
void mp_kji(double** A, double** B, double** C, int n);


#endif
