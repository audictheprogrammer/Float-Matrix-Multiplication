#include "matrix.h"
#include "matrix_mul.h"
#include <math.h>

void benchmark_ijk(double** A, double** B, int n){
    FILE* f = fopen("data/benchmark_order_ijk.txt", "a");
    double** C = zero_matrix(n);

    clock_t initial = clock();
    mp_ijk(A, B, C, n);
    clock_t final = clock();

    double time = ((double) (final - initial)) / CLOCKS_PER_SEC;
    printf("IJK: n = %d time = %f \n", n, time);
    fprintf(f, "%d %f \n", n, time);

    fclose(f);
}

void benchmark_kij(double** A, double** B, int n){
    FILE* f = fopen("data/benchmark_order_kij.txt", "a");
    double** C = zero_matrix(n);

    clock_t initial = clock();
    mp_kij(A, B, C, n);
    clock_t final = clock();

    double time = ((double) (final - initial)) / CLOCKS_PER_SEC;
    printf("KIJ: n = %d time = %f \n", n, time);
    fprintf(f, "%d %f \n", n, time);

    fclose(f);
}

void benchmark_jki(double** A, double** B, int n){
    FILE* f = fopen("data/benchmark_order_jki.txt", "a");
    double** C = zero_matrix(n);

    clock_t initial = clock();
    mp_jki(A, B, C, n);
    clock_t final = clock();

    double time = ((double) (final - initial)) / CLOCKS_PER_SEC;
    printf("JKI: n = %d time = %f \n", n, time);
    fprintf(f, "%d %f \n", n, time);

    fclose(f);
}

void benchmark_ikj(double** A, double** B, int n){
    FILE* f = fopen("data/benchmark_order_ikj.txt", "a");
    double** C = zero_matrix(n);

    clock_t initial = clock();
    mp_ikj(A, B, C, n);
    clock_t final = clock();

    double time = ((double) (final - initial)) / CLOCKS_PER_SEC;
    printf("IKJ: n = %d time = %f \n", n, time);
    fprintf(f, "%d %f \n", n, time);

    fclose(f);
}

void benchmark_jik(double** A, double** B, int n){
    FILE* f = fopen("data/benchmark_order_jik.txt", "a");
    double** C = zero_matrix(n);

    clock_t initial = clock();
    mp_jik(A, B, C, n);
    clock_t final = clock();

    double time = ((double) (final - initial)) / CLOCKS_PER_SEC;
    printf("JIK: n = %d time = %f \n", n, time);
    fprintf(f, "%d %f \n", n, time);

    fclose(f);
}

void benchmark_kji(double** A, double** B, int n){
    FILE* f = fopen("data/benchmark_order_kji.txt", "a");
    double** C = zero_matrix(n);

    clock_t initial = clock();
    mp_kji(A, B, C, n);
    clock_t final = clock();

    double time = ((double) (final - initial)) / CLOCKS_PER_SEC;
    printf("KJI: n = %d time = %f \n\n", n, time);
    fprintf(f, "%d %f \n", n, time);

    fclose(f);
}

void benchmark_loops_order(double p){

    for (int i=8; i<11; i++){
        int n = (int) pow(2, i);
        double**A = random_matrix(n, p);
        double**B = random_matrix(n, p);
        benchmark_ijk(A, B, n);
        benchmark_kij(A, B, n);
        benchmark_jki(A, B, n);
        benchmark_ikj(A, B, n);
        benchmark_jik(A, B, n);
        benchmark_kji(A, B, n);
    }
}

// void benchmark_mod_naive(double** A, double** B, int n){
//     FILE* f = fopen("data/benchmark_kji.txt", "a");
//     double** C = zero_matrix(n);
//
//     clock_t initial = clock();
//     mp_kji(A, B, C, n);
//     clock_t final = clock();
//
//     double time = ((double) (final - initial)) / CLOCKS_PER_SEC;
//     printf("KJI: n = %d time = %f \n\n", n, time);
//     fprintf(f, "%d %f \n", n, time);
//
//     fclose(f);
// }

void benchmark_modulos(double p, double u){

    // for (int i=8; i<11; i++){
    //     int n = (int) pow(2, i);
    //     double**A = random_matrix(n, p);
    //     double**B = random_matrix(n, p);
    //     benchmark_mod_naive(A, B, n);
    //     benchmark_mod_SIMD1(A, B, n);
    //     benchmark_mod_SIMD2(A, B, n);
    //     benchmark_mod_SIMD3(A, B, n);
    //     benchmark_mod_barrett(A, B, n);
    //
    // }
}




int main(){
    // Initialization
    srand(time(NULL));
    double p = pow(2, 26) - 5;
    double u = 1.0 / p;  // u = inv(p)
    u_int32_t u_b = (int) (pow(2, 36) / p);

    // // Testing loops order
    // benchmark_loops_order(p);

    // Testing different modulo
    // benchmark_modulos(p, u);


    return 0;
}
