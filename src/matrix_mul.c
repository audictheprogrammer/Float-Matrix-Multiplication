#include "matrix_mul.h"


double modulo_naive(double a, double p){
    return (double) ((long) a % (long) p);
}

double modulo_SIMD1(double a, double p, double u){
    /* Function 3.1 from SIMD article for floats.
    Hypothesis: Rounding mode = nearest and p < 2^26.
    */
    double b = a * u;
    double c = (double)(int)b;
    double d = a - c * p;
    if (d >= p) return d-p;
    if (d < 0) return d+p;
    return d;
}

double modulo_SIMD2(double a, double p, double u){
    /* Function 3.1 from SIMD article for floats.
    Hypothesis: Rounding mode = up and p < 2^26.
    */
    double b = a * u;
    double c = (double)(int)b;
    double d = a - c * p;
    if (d < 0) return d+p;
    return d;
}

double modulo_SIMD3(double a, double p, double u){
    /* Function 3.1 from SIMD article for floats.
    Hypothesis: Rounding mode = up and p < 2^26.
    */
    double b = a * u;
    double c = (double)(int)b;
    double d = a - c * p;
    // return d + (-(d<0) & (int64_t) p);
    return d + p * (d<0);
}

u_int32_t modulo_Barrett(u_int64_t a, u_int32_t p, u_int32_t u, u_int32_t s, u_int32_t t){
    /* Barrett's modular function for integers.
    Hypothesis: 0 <= a < 2^{32+bitsize_p}, and no assumption on p.
    Returns a % p
    */

    u_int64_t b = a >> s;
    u_int64_t c = (b * u) >> t; // u = 2^(s+t) / p

    u_int32_t res = a - c * p;

    if (res >= p) return res - p;
    return res;
}


int get_bitsize(double p){
    /* Get the bitsize of p.
    Bitsize < 26 because we work with doubles and product.
    */
    double MAX = pow(2, 26);

    for (int i = 26; i > 0; i--){
        if (p > MAX){
            return i+1;  // res = 2^i forall i
        }
        MAX /= 2;
    }
    return 0;
}

int get_blocksize(int b, int n){
    // b: bitsize of p
    // n: size of matrix
    int max_bitsize_double = 53;
    int res = (int) pow(2, max_bitsize_double - 2*b);
    if (res > n){
        return n;
    }
    return res;
}

void mp_naive(double* A, double* B, double* C, int n, double p){
    // Assert C is a zero matrix.
    for (int k=0; k<n; k++){
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                double temp = modulo_naive(A[i*n + k] * B[k*n + j], p);
                C[i*n + j] += temp;
            }
        }
    }

    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            C[i*n + j] = modulo_naive(C[i*n + j], p);
        }
    }

}

void mp_SIMD1(double* A, double* B, double* C, int n, double p, double u){
    // Assert C is a zero matrix
    for (int k=0; k<n; k++){
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                double temp = modulo_SIMD1(A[i*n + k] * B[k*n + j], p, u);
                C[i*n + j] += temp;
            }
        }
    }

    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            C[i*n + j] = modulo_SIMD1(C[i*n + j], p, u);
        }
    }

}

void mp_SIMD2(double* A, double* B, double* C, int n, double p, double u){
    // Assert C is a zero matrix
    for (int k=0; k<n; k++){
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                double temp = modulo_SIMD2(A[i*n + k] * B[k*n + j], p, u);
                C[i*n + j] += temp;
            }
        }
    }

    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            C[i*n + j] = modulo_SIMD2(C[i*n + j], p, u);
        }
    }

}

void mp_SIMD3(double* A, double* B, double* C, int n, double p, double u){
    // Assert C is a zero matrix
    for (int k=0; k<n; k++){
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                double temp = modulo_SIMD3(A[i*n + k] * B[k*n + j], p, u);
                C[i*n + j] += temp;
            }
        }
    }

    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            C[i*n + j] = modulo_SIMD3(C[i*n + j], p, u);
        }
    }

}


void mp_Barrett(double* A, double* B, double* C, int n, double p, u_int32_t u, u_int32_t s, u_int32_t t){
    // Assert C is a zero matrix
    for (int k=0; k<n; k++){
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                double temp = modulo_Barrett(A[i*n + k] * B[k*n + j], p, u, s, t);
                C[i*n + j] += temp;
            }
        }
    }

    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            C[i*n + j] = modulo_Barrett(C[i*n + j], p, u, s, t);
        }
    }

}


// OpenMP
void mp_naive_MP(double* A, double* B, double* C, int n, double p){
    // Assert C is a zero matrix.
    for (int k=0; k<n; k++){
        for (int i=0; i<n; i++){
            #pragma omp parallel for
            for (int j=0; j<n; j++){
                    C[i*n + j] = C[i*n + j] + modulo_naive(A[i*n + k] * B[k*n + j], p);
            }
        }
    }

    for (int i=0; i<n; i++){
        #pragma omp parallel for
        for (int j=0; j<n; j++){
            C[i*n + j] = modulo_naive(C[i*n + j], p);
        }
    }

}

void mp_SIMD1_MP(double* A, double* B, double* C, int n, double p, double u){
    // Assert C is a zero matrix
    for (int k=0; k<n; k++){
        for (int i=0; i<n; i++){
            #pragma omp parallel for
            for (int j=0; j<n; j++){
                C[i*n + j] = C[i*n + j] + modulo_SIMD1(A[i*n + k] * B[k*n + j], p, u);
            }
        }
    }

    for (int i=0; i<n; i++){
        #pragma omp parallel for
        for (int j=0; j<n; j++){
                C[i*n + j] = modulo_SIMD1(C[i*n + j], p, u);
        }
    }

}

void mp_SIMD2_MP(double* A, double* B, double* C, int n, double p, double u){
    // Assert C is a zero matrix
    for (int k=0; k<n; k++){
        for (int i=0; i<n; i++){
            #pragma omp parallel for
            for (int j=0; j<n; j++){
                    C[i*n + j] = C[i*n + j] + modulo_SIMD2(A[i*n + k] * B[k*n + j], p, u);
            }
        }
    }

    for (int i=0; i<n; i++){
        #pragma omp parallel for
        for (int j=0; j<n; j++){
                C[i*n + j] = modulo_SIMD2(C[i*n + j], p, u);
        }
    }

}

void mp_SIMD3_MP(double* A, double* B, double* C, int n, double p, double u){
    // Assert C is a zero matrix
    for (int k=0; k<n; k++){
        for (int i=0; i<n; i++){
            #pragma omp parallel for
            for (int j=0; j<n; j++){
                    C[i*n + j] = C[i*n + j] + modulo_SIMD3(A[i*n + k] * B[k*n + j], p, u);
            }
        }
    }

    for (int i=0; i<n; i++){
        #pragma omp parallel for
        for (int j=0; j<n; j++){
            C[i*n + j] = modulo_SIMD3(C[i*n + j], p, u);
        }
    }

}


void mp_Barrett_MP(double* A, double* B, double* C, int n, double p, u_int32_t u, u_int32_t s, u_int32_t t){
    // Assert C is a zero matrix
    for (int k=0; k<n; k++){
        for (int i=0; i<n; i++){
            #pragma omp parallel for
            for (int j=0; j<n; j++){
                C[i*n + j] = C[i*n + j] + modulo_Barrett(A[i*n + k] * B[k*n + j], p, u, s, t);
            }
        }
    }

    for (int i=0; i<n; i++){
        #pragma omp parallel for
        for (int j=0; j<n; j++){
            C[i*n + j] = modulo_Barrett(C[i*n + j], p, u, s, t);
        }
    }

}



void mp_block(double* A, double* B, double* C, int n, double p, double u, int b){
    /* Compute the product of two matrices using basic block product.
    It allows us to reduce the amount of modulo needed.
    */
    for (int k=0; k<n; k+=b){

        for (int kk=k; kk<k+b; kk++){
            for(int ii=0; ii<n; ii++){
                for (int jj=0; jj<n; jj++){
                    C[ii*n + jj] += A[ii*n + kk] * B[kk*n + jj];
                }
            }
        }

        for (int i=0; i<n*n; i++){
            C[i] = modulo_SIMD3(C[i], p, u);
        }

    }

}

void mp_block_BLAS(double* A, double* B, double* C, int n, double p, double u, int b){
    /* Compute the product of two matrices using OpenBLAS's block product.
    It allows us to reduce the amount of modulo needed.
    */
    for (int k=0; k<n; k+=b){
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    n, n, b, 1, A + k, n, B + n*k, n,
                    1, C, n);


        for (int i=0; i<n*n; i++){
            C[i] = modulo_SIMD3(C[i], p, u);
        }

    }

}

void mp_block_BLAS_MP(double* A, double* B, double* C, int n, double p, double u, int b){
    /* Compute the product of two matrices using OpenBLAS's block product.
    It allows us to reduce the amount of modulo needed.
    */
    for (int k=0; k<n; k+=b){
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    n, n, b, 1, A + k, n, B + n*k, n,
                    1, C, n);

        #pragma omp parallel for
        for (int i=0; i<n*n; i++){
            C[i] = modulo_SIMD3(C[i], p, u);
        }

    }

}




// Comparing loop order. IKJ wins.

// Loop 1
void mp_ijk(double* A, double* B, double* C, int n){
    // Assert C is a zero matrix
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            for (int k=0; k<n; k++){
                C[i*n + j] += A[i*n + k] * B[k*n + j];
            }
        }
    }
}

// Loop 2
void mp_kij(double* A, double* B, double* C, int n){
    // Assert C is a zero matrix
    for (int k=0; k<n; k++){
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                C[i*n + j] += A[i*n + k] * B[k*n + j];
            }
        }
    }
}


// Loop 3
void mp_jki(double* A, double* B, double* C, int n){
    // Assert C is a zero matrix
    for (int j=0; j<n; j++){
        for (int k=0; k<n; k++){
            for (int i=0; i<n; i++){
                C[i*n + j] += A[i*n + k] * B[k*n + j];
            }
        }
    }
}

// Loop 4
void mp_ikj(double* A, double* B, double* C, int n){
    // Assert C is a zero matrix
    for (int i=0; i<n; i++){
        for (int k=0; k<n; k++){
            for (int j=0; j<n; j++){
                C[i*n + j] += A[i*n + k] * B[k*n + j];
            }
        }
    }
}

// Loop 5
void mp_jik(double* A, double* B, double* C, int n){
    // Assert C is a zero matrix
    for (int j=0; j<n; j++){
        for (int i=0; i<n; i++){
            for (int k=0; k<n; k++){
                C[i*n + j] += A[i*n + k] * B[k*n + j];
            }
        }
    }
}

// Loop 6
void mp_kji(double* A, double* B, double* C, int n){
    // Assert C is a zero matrix
    for (int k=0; k<n; k++){
        for (int j=0; j<n; j++){
            for (int i=0; i<n; i++){
                C[i*n + j] += A[i*n + k] * B[k*n + j];
            }
        }
    }
}


void mp_integer(u_int64_t* A, u_int64_t* B, u_int64_t* C, int n, u_int32_t p, u_int32_t u, u_int32_t s, u_int32_t t){
    /* Matrix product for Integers. */
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            u_int64_t h = 0;
            u_int64_t l = 0;
            for (int k = 0; k < n; k++){
                u_int64_t temp = A[i*n + k] * B[j*n + k];
                h += temp >> 32;
                l += temp - (((temp) >> 32) << 32);
            }
            u_int64_t o = 1;
            u_int64_t h_rem = modulo_Barrett(h, p, u, s, t);
            u_int64_t l_rem = modulo_Barrett(l, p, u, s, t);
            C[i*n + j] = modulo_Barrett(h_rem * (o << 32) + l_rem, p, u, s, t);

        }
    }

}

void mp_float(double* A, double* B, double* C, int n, double p, double u, int b){
    /* Matrix product for Floats using: 1 thread only. */
    for (int k=0; k<n; k+=b){
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    n, n, b, 1, A + k, n, B + n*k, n,
                    1, C, n);

        for (int i=0; i<n*n; i++){
            C[i] = modulo_SIMD3(C[i], p, u);
        }

    }

}
