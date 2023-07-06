#ifndef MATRIX_MUL_H
#define MATRIX_MUL_H

#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include <fenv.h>
#include <math.h>

double modulo_naive(double a, double p);
double modulo_SIMD1(double a, double p, double u);
double modulo_SIMD2(double a, double p, double u);
double modulo_SIMD3(double a, double p, double u);
u_int32_t modulo_barrett(u_int64_t a, u_int32_t p, u_int32_t u_b);

void mp_naive(double** A, double** B, double** C, int n, double p);
void mp_SIMD1(double** A, double** B, double** C, int n, double p, double u);
void mp_SIMD2(double** A, double** B, double** C, int n, double p, double u);
void mp_SIMD3(double** A, double** B, double** C, int n, double p, double u);

void mp_ijk(double** A, double** B, double** C, int n);
void mp_kij(double** A, double** B, double** C, int n);
void mp_jki(double** A, double** B, double** C, int n);
void mp_ikj(double** A, double** B, double** C, int n);
void mp_jik(double** A, double** B, double** C, int n);
void mp_kji(double** A, double** B, double** C, int n);


#endif
