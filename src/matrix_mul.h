#ifndef MATRIX_MUL_H
#define MATRIX_MUL_H

#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include <fenv.h>
#include <math.h>
#include <cblas.h>
#include <omp.h>

double modulo_naive(double a, double p);
double modulo_SIMD1(double a, double p, double u);
double modulo_SIMD2(double a, double p, double u);
double modulo_SIMD3(double a, double p, double u);
u_int32_t modulo_Barrett(u_int64_t a, u_int32_t p, u_int32_t u_b);
int get_bitsize(double p);
int get_blocksize(int b, int n);

void mp_naive(double* A, double* B, double* C, int n, double p);
void mp_SIMD1(double* A, double* B, double* C, int n, double p, double u);
void mp_SIMD2(double* A, double* B, double* C, int n, double p, double u);
void mp_SIMD3(double* A, double* B, double* C, int n, double p, double u);
void mp_Barrett(double* A, double* B, double* C, int n, double p, u_int32_t u);

void mp_naive_MP(double* A, double* B, double* C, int n, double p);
void mp_SIMD1_MP(double* A, double* B, double* C, int n, double p, double u);
void mp_SIMD2_MP(double* A, double* B, double* C, int n, double p, double u);
void mp_SIMD3_MP(double* A, double* B, double* C, int n, double p, double u);
void mp_Barrett_MP(double* A, double* B, double* C, int n, double p, u_int32_t u);


void mp_block(double* A, double* B, double* C, int n, double p, double u, int blocksize);
void mp_block_BLAS(double* A, double* B, double* C, int n, double p, double u, int b);


void mp_ijk(double* A, double* B, double* C, int n);
void mp_kij(double* A, double* B, double* C, int n);
void mp_jki(double* A, double* B, double* C, int n);
void mp_ikj(double* A, double* B, double* C, int n);
void mp_jik(double* A, double* B, double* C, int n);
void mp_kji(double* A, double* B, double* C, int n);


#endif
