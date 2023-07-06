#include <stdio.h>
#include <stdlib.h>


double modulo_naive(double a, double p){

    /* Return a % p
    a mod p = frac(a) + (int(a) % p)
    */
    int int_part = ((int) a) % ((int) p);
    double frac_part = a - (double)(int) a;
    printf("Integer part of %f = %d \n", a, int_part);
    printf("Fractional part of %f = %f \n", a, frac_part);
    return frac_part + int_part;

}

double modulo_SIMD1(double a, double p, double u){
    // Function 3.1 of SIMD article
    // double u = 1.0 / p;
    double b = a * u;
    // double c = b;
    int c = (int) b;

    double d = a - c * p;
    if (d >= p) return d-p;
    if (d < 0) return d+p;

    // printf("u = %f\n", u);
    // printf("b = %f\n", b);
    // printf("c = %d\n", c);
    // printf("d = %f\n", d);
    // printf("c*p = %d * %d = %d\n", c, p, c*p);
    return d;
}


void mp_naive(double** A, double** B, double** C, int n, double p){
    // Assert C is a zero matrix
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            for (int k=0; k<n; k++){
                double temp = modulo_naive(A[i][k] * B[k][j], p);
                C[i][j] = modulo_naive(C[i][j] + temp, p);
            }
        }
    }
}

void mp_SIMD1(double** A, double** B, double** C, int n, double p, double u){
    // Assert C is a zero matrix
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            for (int k=0; k<n; k++){
                double temp = modulo_SIMD1(A[i][k] * B[k][j], p, u);
                C[i][j] = modulo_SIMD1(C[i][j] + temp, p, u);
                // 2 modulos
            }
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
