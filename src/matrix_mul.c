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
    // if (d >= p) return d-p;
    if (d >= p) {
        printf("SIMD1: if1 \n");
        return d-p;
    }
    // printf("SIMD1: if1 \n");
    if (d < 0){
        printf("SIMD2: if2 \n");
        return d+p;
    }

    // if (d < 0) return d+p;
    return d;
}

double modulo_SIMD2(double a, double p, double u){
    /* Function 3.1 from SIMD article for floats.
    Hypothesis: Rounding mode = up and p < 2^26.
    */
    double b = a * u;
    double c = (double)(int)b;
    double d = a - c * p;
    // if (d < 0) return d+p;
    if (d < 0){
        printf("d < 0 \n");
        return d+p;
    }
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

u_int32_t modulo_Barrett(u_int64_t a, u_int32_t p, u_int32_t u){
    /* Barrett's modular function for integers.
    Hypothesis: nonnegative numbers, and no assumption on p.
    Returns a % p
    */
    u_int32_t s = 23;
    u_int32_t t = 33;

    u_int64_t b = a >> s;
    u_int64_t c = (b * u) >> t; // u = 2^(s+t) / p

    u_int32_t res = a - c * p;

    if (res >= p) return res - p;
    return res;
}

void mp_naive(double** A, double** B, double** C, int n, double p){
    // Assert C is a zero matrix.
    for (int i=0; i<n; i++){
        for (int k=0; k<n; k++){
            for (int j=0; j<n; j++){
                double temp = modulo_naive(A[i][k] * B[k][j], p);
                C[i][j] += temp;
            }
        }
    }

    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            C[i][j] = modulo_naive(C[i][j], p);
        }
    }

}

void mp_SIMD1(double** A, double** B, double** C, int n, double p, double u){
    // Assert C is a zero matrix
    for (int i=0; i<n; i++){
        for (int k=0; k<n; k++){
            for (int j=0; j<n; j++){
                double temp = modulo_SIMD1(A[i][k] * B[k][j], p, u);
                C[i][j] += temp;
            }
        }
    }

    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            C[i][j] = modulo_SIMD1(C[i][j], p, u);
        }
    }

}

void mp_SIMD2(double** A, double** B, double** C, int n, double p, double u){
    // Assert C is a zero matrix
    for (int i=0; i<n; i++){
        for (int k=0; k<n; k++){
            for (int j=0; j<n; j++){
                double temp = modulo_SIMD2(A[i][k] * B[k][j], p, u);
                C[i][j] += temp;
            }
        }
    }

    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            C[i][j] = modulo_SIMD2(C[i][j], p, u);
        }
    }

}

void mp_SIMD3(double** A, double** B, double** C, int n, double p, double u){
    // Assert C is a zero matrix
    for (int i=0; i<n; i++){
        for (int k=0; k<n; k++){
            for (int j=0; j<n; j++){
                double temp = modulo_SIMD3(A[i][k] * B[k][j], p, u);
                C[i][j] += temp;
            }
        }
    }

    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            C[i][j] = modulo_SIMD3(C[i][j], p, u);
        }
    }

}


void mp_Barrett(double** A, double** B, double** C, int n, double p, u_int32_t u){
    // Assert C is a zero matrix
    for (int i=0; i<n; i++){
        for (int k=0; k<n; k++){
            for (int j=0; j<n; j++){
                double temp = modulo_Barrett(A[i][k] * B[k][j], p, u);
                C[i][j] += temp;
            }
        }
    }

    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            C[i][j] = modulo_Barrett(C[i][j], p, u);
        }
    }

}
//
// void mp_block_2D(double** A, double** B, double** C, int n, double p, double u, int b){
//     /* Compute the product of two 2D matrices using block product.
//     Without OpenBLAS.
//     */
//     for (int i=0; i<n; i+=b){
//         for (int j=0; j<n; j+=b){
//             for (int k=0; k<n; k+=b){
//
//                 for (int ii=i; ii<i+b; ii++){
//                     for (int jj=j; jj<j+b; jj++){
//                         for (int kk=k; kk<k+b; kk++){
//                             C[ii][jj] += A[ii][kk] + B[kk][jj];
//                         }
//                     }
//                 }
//
//             }
//         }
//     }
//
// }

void mp_block(double* A, double* B, double* C, int n, double p, double u, int b){
    /* Compute the product of two 1D matrices using block product.
    Without OpenBLAS.
    */
    for (int i=0; i<n; i+=b){
        for (int j=0; j<n; j+=b){
            for (int k=0; k<n; k+=b){

                for (int ii=i; ii<i+b; ii++){
                    for (int jj=j; jj<j+b; jj++){
                        for (int kk=k; kk<k+b; kk++){
                            C[ii*n + jj] += A[ii*n + kk] * B[kk*n + jj];
                        }
                        C[ii*n + jj] = modulo_SIMD2(C[ii*n + jj], p, u);
                    }
                }
            }
        }
    }

}

void mp_block_BLAS(double* A, double* B, double* C, int n, double p, double u, int b){
    for (int i=0; i<n; i+=b){
        for (int j=0; j<n; j+=b){
            for (int k=0; k<n; k+=b){
                // Block
                // printf("i*n+k = %d \n", i*n+k);
                // printf("k*n+j = %d \n", k*n+j);
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                            b, b, b, 1, A + (i*n + k), n, B + (k*n + j),
                            n, 1, C + (i*n + j), n);
                for (int ii=i; ii<i+b; ii++){
                    for (int jj=j; jj<j+b; jj++){
                        C[ii*n + jj] = modulo_SIMD2(C[ii*n + jj], p, u);
                    }
                }
            }
        }
    }
}


int get_bitsize(double p){
    /* Get the bitsize of p.
    Bitsize < 26 because we work with doubles and product.
    */
    double MAX = pow(2, 26);
    if (p >= MAX)
        return -1;
    if (p <= 0)
        return -1;
    int n = 1;
    for (int i = 1; i<26; i++){
        n = n*2;
        if (p < n)
            return i;
    }
    return 26;
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


// Comparing loop order. IKJ wins.

// Loop 1
void mp_ijk(double** A, double** B, double** C, int n){
    // Assert C is a zero matrix
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            for (int k=0; k<n; k++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}


// Loop 2
void mp_kij(double** A, double** B, double** C, int n){
    // Assert C is a zero matrix
    for (int k=0; k<n; k++){
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}


// Loop 3
void mp_jki(double** A, double** B, double** C, int n){
    // Assert C is a zero matrix
    for (int j=0; j<n; j++){
        for (int k=0; k<n; k++){
            for (int i=0; i<n; i++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// Loop 4
void mp_ikj(double** A, double** B, double** C, int n){
    // Assert C is a zero matrix
    for (int i=0; i<n; i++){
        for (int k=0; k<n; k++){
            for (int j=0; j<n; j++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// Loop 5
void mp_jik(double** A, double** B, double** C, int n){
    // Assert C is a zero matrix
    for (int j=0; j<n; j++){
        for (int i=0; i<n; i++){
            for (int k=0; k<n; k++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// Loop 6
void mp_kji(double** A, double** B, double** C, int n){
    // Assert C is a zero matrix
    for (int k=0; k<n; k++){
        for (int j=0; j<n; j++){
            for (int i=0; i<n; i++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}
