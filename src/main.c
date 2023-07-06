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

void benchmark_mod_naive(double** A, double** B, int n, double p){
    FILE* f = fopen("data/benchmark_modulo_naive.txt", "a");
    double** C = zero_matrix(n);

    clock_t initial = clock();
    mp_naive(A, B, C, n, p);
    clock_t final = clock();

    double time = ((double) (final - initial)) / CLOCKS_PER_SEC;
    printf("Modulo Naive: n = %d time = %f \n\n", n, time);
    fprintf(f, "%d %f \n", n, time);

    fclose(f);
}

void benchmark_mod_SIMD1(double** A, double** B, int n, double p, double u){
    FILE* f = fopen("data/benchmark_modulo_SIMD1.txt", "a");
    double** C = zero_matrix(n);

    clock_t initial = clock();
    mp_SIMD1(A, B, C, n, p, u);
    clock_t final = clock();

    double time = ((double) (final - initial)) / CLOCKS_PER_SEC;
    printf("Modulo SIMD1: n = %d time = %f \n\n", n, time);
    fprintf(f, "%d %f \n", n, time);

    fclose(f);
}

void benchmark_mod_SIMD2(double** A, double** B, int n, double p, double u){
    FILE* f = fopen("data/benchmark_modulo_SIMD2.txt", "a");
    double** C = zero_matrix(n);

    clock_t initial = clock();
    mp_SIMD2(A, B, C, n, p, u);
    clock_t final = clock();

    double time = ((double) (final - initial)) / CLOCKS_PER_SEC;
    printf("Modulo SIMD2: n = %d time = %f \n\n", n, time);
    fprintf(f, "%d %f \n", n, time);

    fclose(f);
}

void benchmark_mod_SIMD3(double** A, double** B, int n, double p, double u){
    FILE* f = fopen("data/benchmark_modulo_SIMD3.txt", "a");
    double** C = zero_matrix(n);

    clock_t initial = clock();
    mp_SIMD3(A, B, C, n, p, u);
    clock_t final = clock();

    double time = ((double) (final - initial)) / CLOCKS_PER_SEC;
    printf("Modulo SIMD3: n = %d time = %f \n\n", n, time);
    fprintf(f, "%d %f \n", n, time);

    fclose(f);
}

void benchmark_modulos(double p, double u){

    for (int i=8; i<11; i++){
        int n = (int) pow(2, i);
        double**A = random_matrix(n, p);
        double**B = random_matrix(n, p);
        benchmark_mod_naive(A, B, n, p);  // Worst
        benchmark_ikj(A, B, n);  // Best
        benchmark_mod_SIMD1(A, B, n, p, u);
        benchmark_mod_SIMD2(A, B, n, p, u);
        benchmark_mod_SIMD3(A, B, n, p, u);

    }

}

void clean_file_loops(){
    char noms[6][64] = {"data/benchmark_order_ijk.txt", "data/benchmark_order_ikj.txt",\
                    "data/benchmark_order_jik.txt", "data/benchmark_order_jki.txt",\
                    "data/benchmark_order_kij.txt", "data/benchmark_order_kji.txt"};

    for (int i=0; i<6; i++){
        FILE* f = fopen(noms[i], "w");
        fclose(f);
    }
}

void clean_file_modulos(){
    char noms[5][64] = {"data/benchmark_modulo_naive.txt", "data/benchmark_order_ijk.txt",\
                    "data/benchmark_modulo_SIMD1.txt", "data/benchmark_modulo_SIMD2.txt",\
                    "data/benchmark_modulo_SIMD3.txt"};

    for (int i=0; i<5; i++){
        FILE* f = fopen(noms[i], "w");
        fclose(f);
    }
}


int main(){
    // Initialization
    srand(time(NULL));
    double p = pow(2, 26) - 5;
    double u = 1.0 / p;  // u = inv(p)
    u_int32_t u_b = (int) (pow(2, 36) / p);


    // // // Testing loops order
    // clean_file_loops();
    // benchmark_loops_order(p);

    // Testing different modulo
    clean_file_modulos();
    benchmark_modulos(p, u);


    return 0;
}
