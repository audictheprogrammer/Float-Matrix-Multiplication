#include "matrix.h"
#include "matrix_mul.h"
#include <math.h>

void benchmark_ijk(char* filename, float** A, float** B, int n){
    FILE* f = fopen(filename, "a");
    float** C = zero_matrix(n);

    clock_t initial = clock();
    mp_ijk(A, B, C, n);
    clock_t final = clock();

    double time = ((double) (final - initial)) / CLOCKS_PER_SEC;
    printf("IJK: n = %d time = %f \n", n, time);
    fprintf(f, "IJK: n = %d time = %f \n", n, time);

    fclose(f);
}

void benchmark_kij(char* filename, float** A, float** B, int n){
    FILE* f = fopen(filename, "a");
    float** C = zero_matrix(n);

    clock_t initial = clock();
    mp_kij(A, B, C, n);
    clock_t final = clock();

    double time = ((double) (final - initial)) / CLOCKS_PER_SEC;
    printf("KIJ: n = %d time = %f \n", n, time);
    fprintf(f, "KIJ: n = %d time = %f \n", n, time);

    fclose(f);
}

void benchmark_jki(char* filename, float** A, float** B, int n){
    FILE* f = fopen(filename, "a");
    float** C = zero_matrix(n);

    clock_t initial = clock();
    mp_jki(A, B, C, n);
    clock_t final = clock();

    double time = ((double) (final - initial)) / CLOCKS_PER_SEC;
    printf("JKI: n = %d time = %f \n", n, time);
    fprintf(f, "JKI: n = %d time = %f \n", n, time);

    fclose(f);
}

void benchmark_ikj(char* filename, float** A, float** B, int n){
    FILE* f = fopen(filename, "a");
    float** C = zero_matrix(n);

    clock_t initial = clock();
    mp_ikj(A, B, C, n);
    clock_t final = clock();

    double time = ((double) (final - initial)) / CLOCKS_PER_SEC;
    printf("IKJ: n = %d time = %f \n", n, time);
    fprintf(f, "IKJ: n = %d time = %f \n", n, time);

    fclose(f);
}

void benchmark_jik(char* filename, float** A, float** B, int n){
    FILE* f = fopen(filename, "a");
    float** C = zero_matrix(n);

    clock_t initial = clock();
    mp_jik(A, B, C, n);
    clock_t final = clock();

    double time = ((double) (final - initial)) / CLOCKS_PER_SEC;
    printf("JIK: n = %d time = %f \n", n, time);
    fprintf(f, "JIK: n = %d time = %f \n", n, time);

    fclose(f);
}

void benchmark_kji(char* filename, float** A, float** B, int n){
    FILE* f = fopen(filename, "a");
    float** C = zero_matrix(n);

    clock_t initial = clock();
    mp_kji(A, B, C, n);
    clock_t final = clock();

    double time = ((double) (final - initial)) / CLOCKS_PER_SEC;
    printf("KJI: n = %d time = %f \n\n", n, time);
    fprintf(f, "KJI: n = %d time = %f \n\n", n, time);

    fclose(f);
}

void benchmark_loops_order(char* name){
    char* filename = strdup(name);
    int p = (int) pow(2, 26)- 5;  // A prime number

    for (int i=8; i<11; i++){
        int n = (int) pow(2, i);
        float**A = random_matrix(n, p);
        float**B = random_matrix(n, p);
        float**C = zero_matrix(n);
        benchmark_ijk(filename, A, B, n);
        benchmark_kij(filename, A, B, n);
        benchmark_jki(filename, A, B, n);
        benchmark_ikj(filename, A, B, n);
        benchmark_jik(filename, A, B, n);
        benchmark_kji(filename, A, B, n);
    }
    free(filename);
}



int main(){
    // Testing loops order
    // benchmark_loops_order("benchmark_loops_order.txt");

    int p = (int) pow(2, 26) - 5;
    int n = 2;
    // printf("Single test\n");
    // float**A = random_matrix(n, p);
    // float**B = random_matrix(n, p);
    // float**C = zero_matrix(n);
    // mp_ijk(A, B, C, n);
    //
    // print_matrix(A, n);
    // print_matrix(B, n);
    // print_matrix(C, n);

    float a = 1024.2;
    printf("%f mod %d = %f \n", a, p, modulo_naive(a, p));


    return 0;
}
