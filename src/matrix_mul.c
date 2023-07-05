#include <stdio.h>
#include <stdlib.h>


float modulo_naive(float a, int p){

    /* Return a % p
    a mod p = frac(a) + (int(a) % p)
    */
    int int_part = ((int) a) % p;
    float frac_part = a - (float)(int) a;
    printf("Integer part of %f = %d \n", a, int_part);
    printf("Fractional part of %f = %f \n", a, frac_part);
    return frac_part + int_part;

}

float modulo_SIMD1(float a, int p){
    // Function 3.1 of SIMD article
    float u = 1.0 / p;
    float b = a * u;
    int c = (int) b;

    float d = a - c * p;
    if (d >= p) return d-p;
    if (d < 0) return d+p;

    printf("u = %f\n", u);
    printf("b = %f\n", b);
    printf("c = %d\n", c);
    printf("d = %f\n", d);
    printf("c*p = %d * %d = %d\n", c, p, c*p);
    return d;
}


void mp_naive(float** A, float** B, float** C, int n, int p){
    // Assert C is a zero matrix
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            for (int k=0; k<n; k++){
                C[i][j] += A[i][k] * B[k][j];
                C[i][j] = modulo_naive(C[i][j], p);
            }
        }
    }
}


// Comparing loop order. KIJ wins.

// Loop 1
void mp_ijk(float** A, float** B, float** C, int n){
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
void mp_kij(float** A, float** B, float** C, int n){
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
void mp_jki(float** A, float** B, float** C, int n){
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
void mp_ikj(float** A, float** B, float** C, int n){
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
void mp_jik(float** A, float** B, float** C, int n){
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
void mp_kji(float** A, float** B, float** C, int n){
    // Assert C is a zero matrix
    for (int k=0; k<n; k++){
        for (int j=0; j<n; j++){
            for (int i=0; i<n; i++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}
