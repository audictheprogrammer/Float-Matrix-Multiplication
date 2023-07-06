#include "matrix_mul.h"


double modulo_naive(double a, double p){
    return (double) ((long) a % (long) p);
}

double modulo_SIMD1(double a, double p, double u){
    /* Function 3.1 of SIMD article
    Hypothesis: Rounding mode = down and p < 2^26.
    */
    double b = a * u;
    double c = (double)(int)b;
    double d = a - c * p;
    if (d >= p) return d-p;
    if (d < 0) return d+p;

    // printf("u = %f\n", u);
    // printf("b = %f\n", b);
    // printf("c = %f\n", c);
    // printf("d = %f\n", d);
    // printf("c*p = %d * %d = %d\n", c, p, c*p);
    return d;
}

double modulo_SIMD2(double a, double p, double u){
    /* Function 3.1 of SIMD article
    Hypothesis: Rounding mode = up and p < 2^26.
    */
    double b = a * u;
    double c = rint(b);  //  The lib fenv.h has function double rint(double d)
    double d = a - c * p;
    if (d < 0) return d+p;
    return d;
}

double modulo_SIMD3(double a, double p, double u){
    /* Function 3.1 of SIMD article
    Hypothesis: Rounding mode = up and p < 2^26.
    */
    double b = a * u;
    double c = rint(b);
    // The lib fenv.h has a function
    // double rint(double d)
    double d = a - c * p;
    if (d < 0) return d+p;
    int cond1 = d<0;
    // return cond1*(d+p) + (!cond1)*d;
    // return cond1&&(d+p) || (!cond1)&&d;
    return d + (-(d<0) && p);
}

u_int32_t modulo_barrett(u_int64_t a, u_int32_t p, u_int32_t u_b){
    /* Function 3.1 of SIMD article
    Hypothesis: nonnegative numbers, and no assumption on p.
    Returns a % p
    */
    u_int32_t s = 23;
    u_int32_t t = 33;
    // u_int32_t s = 30;
    // u_int32_t t = 32;

    u_int64_t b = a >> s;
    u_int64_t c = (b * u_b) >> t; // 4294967284 = (2^62 / p)

    u_int32_t res = a - c * p;
    while (res >= p){
        res -= p;
    }
    // (res >= p) ? res -= p : res;
    // (res >= p) ? res -= p : res;
    return res;
}

void mp_naive(double** A, double** B, double** C, int n, double p){
    // Assert C is a zero matrix
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            for (int k=0; k<n; k++){
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
        for (int j=0; j<n; j++){
            for (int k=0; k<n; k++){
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
        for (int j=0; j<n; j++){
            for (int k=0; k<n; k++){
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
        for (int j=0; j<n; j++){
            for (int k=0; k<n; k++){
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

// Comparing loop order. KIJ wins.

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
