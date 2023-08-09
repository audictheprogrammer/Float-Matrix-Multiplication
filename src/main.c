#include "matrix.h"
#include "matrix_mul.h"
#include <math.h>


// Matrix Product with different order of loops.
double benchmark_ijk(double* A, double* B, int n){
    double* C = zero_matrix_1D(n);

    clock_t initial = clock();
    mp_ijk(A, B, C, n);
    clock_t final = clock();

    delete_matrix_1D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
}

double benchmark_kij(double* A, double* B, int n){
    double* C = zero_matrix_1D(n);

    clock_t initial = clock();
    mp_kij(A, B, C, n);
    clock_t final = clock();

    delete_matrix_1D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
}

double benchmark_jki(double* A, double* B, int n){
    double* C = zero_matrix_1D(n);

    clock_t initial = clock();
    mp_jki(A, B, C, n);
    clock_t final = clock();

    delete_matrix_1D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
}

double benchmark_ikj(double* A, double* B, int n){
    double* C = zero_matrix_1D(n);

    clock_t initial = clock();
    mp_ikj(A, B, C, n);
    clock_t final = clock();

    delete_matrix_1D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
}

double benchmark_jik(double* A, double* B, int n){
    double* C = zero_matrix_1D(n);

    clock_t initial = clock();
    mp_jik(A, B, C, n);
    clock_t final = clock();

    delete_matrix_1D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
}

double benchmark_kji(double* A, double* B, int n){
    double* C = zero_matrix_1D(n);

    clock_t initial = clock();
    mp_kji(A, B, C, n);
    clock_t final = clock();

    delete_matrix_1D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
}


// Matrix Product with modulos.

double benchmark_mod_naive(double* A, double* B, int n, double p){
    double* C = zero_matrix_1D(n);

    clock_t initial = clock();
    mp_naive(A, B, C, n, p);
    clock_t final = clock();

    delete_matrix_1D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
}

double benchmark_mod_SIMD1(double* A, double* B, int n, double p, double u){
    double* C = zero_matrix_1D(n);

    clock_t initial = clock();
    mp_SIMD1(A, B, C, n, p, u);
    clock_t final = clock();

    delete_matrix_1D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
}

double benchmark_mod_SIMD2(double* A, double* B, int n, double p, double u){
    double* C = zero_matrix_1D(n);

    clock_t initial = clock();
    mp_SIMD2(A, B, C, n, p, u);
    clock_t final = clock();

    delete_matrix_1D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
}

 double benchmark_mod_SIMD3(double* A, double* B, int n, double p, double u){
    double* C = zero_matrix_1D(n);

    clock_t initial = clock();
    mp_SIMD3(A, B, C, n, p, u);
    clock_t final = clock();

    delete_matrix_1D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
}

double benchmark_mod_Barrett(double* A, double* B, int n, double p, double u, u_int32_t s, u_int32_t t){
    double* C = zero_matrix_1D(n);

    clock_t initial = clock();
    mp_Barrett(A, B, C, n, p, u, s, t);
    clock_t final = clock();

    delete_matrix_1D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
}


// Matrix product with modulos parallelized.
double benchmark_mod_MP_naive(double* A, double* B, int n, double p){
    double* C = zero_matrix_1D(n);
    struct timespec initial, final;
    double elapsed;

    clock_gettime(CLOCK_MONOTONIC, &initial);
    mp_naive_MP(A, B, C, n, p);
    clock_gettime(CLOCK_MONOTONIC, &final);

    elapsed = (final.tv_sec - initial.tv_sec);
    elapsed += (final.tv_nsec - initial.tv_nsec) / 1000000000.0;

    delete_matrix_1D(&C, n);
    return elapsed;
}

double benchmark_mod_MP_SIMD1(double* A, double* B, int n, double p, double u){
    double* C = zero_matrix_1D(n);
    struct timespec initial, final;
    double elapsed;

    clock_gettime(CLOCK_MONOTONIC, &initial);
    mp_SIMD1_MP(A, B, C, n, p, u);
    clock_gettime(CLOCK_MONOTONIC, &final);

    elapsed = (final.tv_sec - initial.tv_sec);
    elapsed += (final.tv_nsec - initial.tv_nsec) / 1000000000.0;

    delete_matrix_1D(&C, n);
    return elapsed;
}

double benchmark_mod_MP_SIMD2(double* A, double* B, int n, double p, double u){
    double* C = zero_matrix_1D(n);
    struct timespec initial, final;
    double elapsed;

    clock_gettime(CLOCK_MONOTONIC, &initial);
    mp_SIMD2_MP(A, B, C, n, p, u);
    clock_gettime(CLOCK_MONOTONIC, &final);

    elapsed = (final.tv_sec - initial.tv_sec);
    elapsed += (final.tv_nsec - initial.tv_nsec) / 1000000000.0;

    delete_matrix_1D(&C, n);
    return elapsed;
}

 double benchmark_mod_MP_SIMD3(double* A, double* B, int n, double p, double u){
     double* C = zero_matrix_1D(n);
     struct timespec initial, final;
     double elapsed;

     clock_gettime(CLOCK_MONOTONIC, &initial);
     mp_SIMD2_MP(A, B, C, n, p, u);
     clock_gettime(CLOCK_MONOTONIC, &final);

     elapsed = (final.tv_sec - initial.tv_sec);
     elapsed += (final.tv_nsec - initial.tv_nsec) / 1000000000.0;

     delete_matrix_1D(&C, n);
     return elapsed;
}

double benchmark_mod_MP_Barrett(double* A, double* B, int n, double p, double u, u_int32_t s, u_int32_t t){
    double* C = zero_matrix_1D(n);
    struct timespec initial, final;
    double elapsed;

    clock_gettime(CLOCK_MONOTONIC, &initial);
    mp_Barrett_MP(A, B, C, n, p, u, s, t);
    clock_gettime(CLOCK_MONOTONIC, &final);

    elapsed = (final.tv_sec - initial.tv_sec);
    elapsed += (final.tv_nsec - initial.tv_nsec) / 1000000000.0;

    delete_matrix_1D(&C, n);
    return elapsed;
}


// Blocks Product

double benchmark_blocks_NoBLAS(double* A, double* B, int n, double p, double u, int b){
    double* C = zero_matrix_1D(n);

    clock_t initial = clock();
    mp_block(A, B, C, n, p, u, b);
    clock_t final = clock();

    delete_matrix_1D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
}

double benchmark_blocks_BLAS(double* A, double* B, int n, double p, double u, int b){
    double* C = zero_matrix_1D(n);
    struct timespec initial, final;
    double elapsed;

    clock_gettime(CLOCK_MONOTONIC, &initial);
    mp_block_BLAS(A, B, C, n, p, u, b);
    clock_gettime(CLOCK_MONOTONIC, &final);

    elapsed = (final.tv_sec - initial.tv_sec);
    elapsed += (final.tv_nsec - initial.tv_nsec) / 1000000000.0;

    delete_matrix_1D(&C, n);
    return elapsed;
}


double benchmark_blocks_BLAS_MP(double* A, double* B, int n, double p, double u, int b){
    double* C = zero_matrix_1D(n);
    struct timespec initial, final;
    double elapsed;

    clock_gettime(CLOCK_MONOTONIC, &initial);
    mp_block_BLAS_MP(A, B, C, n, p, u, b);
    clock_gettime(CLOCK_MONOTONIC, &final);

    elapsed = (final.tv_sec - initial.tv_sec);
    elapsed += (final.tv_nsec - initial.tv_nsec) / 1000000000.0;

    delete_matrix_1D(&C, n);
    return elapsed;
}

// Comparing float and integer performance
double benchmark_float(double* A, double* B, int n, double p, double u, int b){
    double* C = zero_matrix_1D(n);
    struct timespec initial, final;
    double elapsed;

    clock_gettime(CLOCK_MONOTONIC, &initial);
    mp_float(A, B, C, n, p, u, b);
    clock_gettime(CLOCK_MONOTONIC, &final);

    elapsed = (final.tv_sec - initial.tv_sec);
    elapsed += (final.tv_nsec - initial.tv_nsec) / 1000000000.0;

    delete_matrix_1D(&C, n);
    return elapsed;
}

double benchmark_integer(u_int64_t* A, u_int64_t* B, int n, u_int32_t p, double u, u_int32_t s, u_int32_t t){
    u_int64_t* C = zero_matrix_1D_integer(n);
    struct timespec initial, final;
    double elapsed;

    clock_gettime(CLOCK_MONOTONIC, &initial);
    mp_integer(A, B, C, n, p, u, s, t);
    clock_gettime(CLOCK_MONOTONIC, &final);

    elapsed = (final.tv_sec - initial.tv_sec);
    elapsed += (final.tv_nsec - initial.tv_nsec) / 1000000000.0;

    delete_matrix_1D_integer(&C, n);
    return elapsed;
}



void write_benchmark_time(char* filename, char* text, int n, double time){
    FILE* f = fopen(filename, "a");
    printf("%s: n = %d time = %f \n", text,n, time);
    fprintf(f, "%d %f \n", n, time);
    fclose(f);
}


void benchmark_loops_order(double p){
    /* Benchmarking the order of loops.
    The most efficient one is IKJ.
    */
    int m = 1;  // Executes m times each algo
    for (int i=8; i<12; i++){
        int n = (int) pow(2, i);
        double sum_ijk = 0;
        double sum_kij = 0;
        double sum_jki = 0;
        double sum_ikj = 0;
        double sum_jik = 0;
        double sum_kji = 0;

        for (int j=0; j<m; j++){
            double* A = random_matrix_1D(n, p);
            double* B = random_matrix_1D(n, p);
            sum_ijk += benchmark_ijk(A, B, n);
            sum_kij += benchmark_kij(A, B, n);
            sum_jki += benchmark_jki(A, B, n);
            sum_ikj += benchmark_ikj(A, B, n);
            sum_jik += benchmark_jik(A, B, n);
            sum_kji += benchmark_kji(A, B, n);

            delete_matrix_1D(&A, n);
            delete_matrix_1D(&B, n);
        }
        printf("\n");
        write_benchmark_time("data/benchmark_order_ijk.txt", "IJK", n, sum_ijk/m);
        write_benchmark_time("data/benchmark_order_kij.txt", "KIJ", n, sum_kij/m);
        write_benchmark_time("data/benchmark_order_jki.txt", "JKI", n, sum_jki/m);
        write_benchmark_time("data/benchmark_order_ikj.txt", "IKJ", n, sum_ikj/m);
        write_benchmark_time("data/benchmark_order_jik.txt", "JIK", n, sum_jik/m);
        write_benchmark_time("data/benchmark_order_kji.txt", "KJI", n, sum_kji/m);

    }
}

void benchmark_modulos(double p, double u, double u_overline, double u_b, u_int32_t s, u_int32_t t){
    /* Benchmarking different modulos.
    The 3 SIMD are the most efficient.
    */
    int m = 1;  // Executes m times each algo
    for (int i=8; i<12; i++){
        int n = (int) pow(2, i);
        double sum_mod_naive = 0;
        double sum_mod_SIMD1 = 0;
        double sum_mod_SIMD2 = 0;
        double sum_mod_SIMD3 = 0;
        double sum_mod_Barrett = 0;

        for (int j=0; j<m; j++){
            double* A = random_matrix_1D(n, p);
            double* B = random_matrix_1D(n, p);
            sum_mod_naive += benchmark_mod_naive(A, B, n, p);
            sum_mod_SIMD1 += benchmark_mod_SIMD1(A, B, n, p, u);
            sum_mod_SIMD2 += benchmark_mod_SIMD2(A, B, n, p, u_overline);
            sum_mod_SIMD3 += benchmark_mod_SIMD3(A, B, n, p, u_overline);
            sum_mod_Barrett += benchmark_mod_Barrett(A, B, n, p, u_b, s, t);

            delete_matrix_1D(&A, n);
            delete_matrix_1D(&B, n);
        }

        printf("\n");
        write_benchmark_time("data/benchmark_modulo_naive.txt", "Mod Naive", n, sum_mod_naive/m);
        write_benchmark_time("data/benchmark_modulo_SIMD1.txt", "Mod SIMD1", n, sum_mod_SIMD1/m);
        write_benchmark_time("data/benchmark_modulo_SIMD2.txt", "Mod SIMD2", n, sum_mod_SIMD2/m);
        write_benchmark_time("data/benchmark_modulo_SIMD3.txt", "Mod SIMD3", n, sum_mod_SIMD3/m);
        write_benchmark_time("data/benchmark_modulo_Barrett.txt", "Mod Barrett", n, sum_mod_Barrett/m);

    }
}

void benchmark_modulos_MP(double p, double u, double u_overline, double u_b, u_int32_t s, u_int32_t t){
    /* Benchmarking different modulos.
    */
    int m = 1;  // Executes m times each algo
    for (int i=8; i<12; i++){
        int n = (int) pow(2, i);

        double sum_mod_MP_naive = 0;
        double sum_mod_MP_SIMD1 = 0;
        double sum_mod_MP_SIMD2 = 0;
        double sum_mod_MP_SIMD3 = 0;
        double sum_mod_MP_Barrett = 0;

        for (int j=0; j<m; j++){
            double* A = random_matrix_1D(n, p);
            double* B = random_matrix_1D(n, p);
            sum_mod_MP_naive += benchmark_mod_MP_naive(A, B, n, p);
            sum_mod_MP_SIMD1 += benchmark_mod_MP_SIMD1(A, B, n, p, u);
            sum_mod_MP_SIMD2 += benchmark_mod_MP_SIMD2(A, B, n, p, u_overline);
            sum_mod_MP_SIMD3 += benchmark_mod_MP_SIMD3(A, B, n, p, u_overline);
            sum_mod_MP_Barrett += benchmark_mod_MP_Barrett(A, B, n, p, u_b, s, t);

            delete_matrix_1D(&A, n);
            delete_matrix_1D(&B, n);
        }

        printf("\n");
        write_benchmark_time("data/benchmark_modulo_MP_naive.txt", "Mod MP Naive", n, sum_mod_MP_naive/m);
        write_benchmark_time("data/benchmark_modulo_MP_SIMD1.txt", "Mod MP SIMD1", n, sum_mod_MP_SIMD1/m);
        write_benchmark_time("data/benchmark_modulo_MP_SIMD2.txt", "Mod MP SIMD2", n, sum_mod_MP_SIMD2/m);
        write_benchmark_time("data/benchmark_modulo_MP_SIMD3.txt", "Mod MP SIMD3", n, sum_mod_MP_SIMD3/m);
        write_benchmark_time("data/benchmark_modulo_MP_Barrett.txt", "Mod MP Barrett", n, sum_mod_MP_Barrett/m);

    }
}


void benchmark_blocks(double p, double u_overline){
    /* Benchmarking different modulos.
    */
    int m = 1;  // Executes m times each algo
    openblas_set_num_threads(1);  // 8 is slower than 1.
    for (int i=8; i<13; i++){
        int n = (int) pow(2, i);
        int b = get_blocksize(get_bitsize(p), n);

        double sum_blocks = 0;
        double sum_blocks_BLAS = 0;

        for (int j=0; j<m; j++){
            double* A = random_matrix_1D(n, p);
            double* B = random_matrix_1D(n, p);
            sum_blocks += benchmark_blocks_NoBLAS(A, B, n, p, u_overline, b);
            sum_blocks_BLAS += benchmark_blocks_BLAS(A, B, n, p, u_overline, b);

            delete_matrix_1D(&A, n);
            delete_matrix_1D(&B, n);
        }

        printf("\n");
        write_benchmark_time("data/benchmark_blocks_NoBLAS.txt", "Block", n, sum_blocks/m);
        write_benchmark_time("data/benchmark_blocks_BLAS.txt", "Block BLAS", n, sum_blocks_BLAS/m);

    }
}


void benchmark_blocks_MP(double p, double u_overline){
    /* Benchmarking different modulos.
    */
    int m = 1;  // Executes m times each algo
    openblas_set_num_threads(1);
    for (int i=8; i<13; i++){
        int n = (int) pow(2, i);
        int b = get_blocksize(get_bitsize(p), n);

        double sum_blocks_BLAS_MP = 0;

        for (int j=0; j<m; j++){
            double* A = random_matrix_1D(n, p);
            double* B = random_matrix_1D(n, p);
            sum_blocks_BLAS_MP += benchmark_blocks_BLAS_MP(A, B, n, p, u_overline, b);

            delete_matrix_1D(&A, n);
            delete_matrix_1D(&B, n);
        }

        printf("\n");
        write_benchmark_time("data/benchmark_blocks_BLAS_MP.txt", "Block BLAS MP", n, sum_blocks_BLAS_MP/m);

    }
}


void benchmark_float_integer(double p, double u_overline, double u_b, u_int32_t s, u_int32_t t){
    /* Benchmarking the integer version of mp.
    */
    int m = 1;  // Executes m times each algo
    openblas_set_num_threads(1);
    for (int i=8; i<13; i++){
        int n = (int) pow(2, i);
        int b = get_blocksize(get_bitsize(p), n);

        double sum_benchmark_float = 0;
        double sum_benchmark_integer = 0;

        for (int j=0; j<m; j++){
            double* A = random_matrix_1D(n, p);
            double* B = random_matrix_1D(n, p);
            u_int64_t* A_int = convert_float_to_integer(A, n);
            u_int64_t* B_int = convert_float_to_integer(B, n);

            sum_benchmark_float += benchmark_float(A, B, n, p, u_overline, b);
            sum_benchmark_integer += benchmark_integer(A_int, B_int, n, p, u_b, s, t);

            delete_matrix_1D(&A, n);
            delete_matrix_1D(&B, n);
            delete_matrix_1D_integer(&A_int, n);
            delete_matrix_1D_integer(&B_int, n);
        }

        printf("\n");

        char filename_float[256] = "\0";
        sprintf(filename_float, "data/benchmark_float_%d.txt", get_bitsize(p));

        char filename_integer[256] = "\0";
        sprintf(filename_integer, "data/benchmark_integer_%d.txt", get_bitsize(p));

        write_benchmark_time(filename_float, "Float mp", n, sum_benchmark_float/m);
        write_benchmark_time(filename_integer, "Integer mp", n, sum_benchmark_integer/m);
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
    char noms[5][64] = {"data/benchmark_modulo_naive.txt", "data/benchmark_modulo_SIMD1.txt",\
     "data/benchmark_modulo_SIMD2.txt", "data/benchmark_modulo_SIMD3.txt",
     "data/benchmark_modulo_Barrett.txt"};

    for (int i=0; i<5; i++){
        FILE* f = fopen(noms[i], "w");
        fclose(f);
    }
}

void clean_file_modulos_MP(){
    char noms[5][64] = {"data/benchmark_modulo_MP_naive.txt", "data/benchmark_modulo_MP_SIMD1.txt",\
     "data/benchmark_modulo_MP_SIMD2.txt", "data/benchmark_modulo_MP_SIMD3.txt",
     "data/benchmark_modulo_MP_Barrett.txt"};

    for (int i=0; i<5; i++){
        FILE* f = fopen(noms[i], "w");
        fclose(f);
    }
}

void clean_file_blocks(){
    char noms[2][64] = {"data/benchmark_blocks_BLAS.txt", "data/benchmark_blocks_NoBLAS.txt"};

    for (int i=0; i<2; i++){
        FILE *f = fopen(noms[i], "w");
        fclose(f);
    }

}

void clean_file_blocks_MP(){
    char noms[1][64] = {"data/benchmark_blocks_BLAS_MP.txt"};

    for (int i=0; i<1; i++){
        FILE *f = fopen(noms[i], "w");
        fclose(f);
    }

}

void clean_file_float_integer(){
    char noms[10][64] = {"data/benchmark_integer_18.txt", "data/benchmark_float_18.txt",
                         "data/benchmark_integer_20.txt", "data/benchmark_float_20.txt",
                         "data/benchmark_integer_22.txt", "data/benchmark_float_22.txt",
                         "data/benchmark_integer_24.txt", "data/benchmark_float_24.txt",
                         "data/benchmark_integer_26.txt", "data/benchmark_float_26.txt"};

    for (int i=0; i<10; i++){
        FILE *f = fopen(noms[i], "w");
        fclose(f);
    }

}


int main(){
    // Initialization
    srand(time(NULL));
    double p = pow(2, 26) - 5;
    // double p = pow(2, 24) - 3;
    // double p = pow(2, 22) - 3;
    // double p = pow(2, 20) - 3;
    // double p = pow(2, 18) - 5;

    // Precomputed constants for Modular functions
    double u = 1.0 / p;  // Constant for SIMD
    fesetround(FE_UPWARD);
    double u_overline = 1.0 / p;  // Constant for SIMD2 and SIMD3
    fesetround(FE_TONEAREST);
    // Constants for Barrett reductionm, which is the integer version of SIMD
    u_int32_t t = 32;
    u_int32_t s = 30 + get_bitsize(p) - t;
    u_int32_t u_b = (int) (pow(2, s+t) / p);


    /* Benchmarking the order of loop. */
    clean_file_loops();
    benchmark_loops_order(p);

    /* Benchmarking different modulo techniques. */
    clean_file_modulos();
    benchmark_modulos(p, u, u_overline, u_b, s, t);

    /* Benchmarking parallelism apporach with previous modulo techniques. */
    clean_file_modulos_MP();
    benchmark_modulos_MP(p, u, u_overline, u_b, s, t);

    /* Benchmarking block product approach. */
    clean_file_blocks();
    benchmark_blocks(p, u_overline);

    /* Benchmarking the final mixed implementation. */
    clean_file_blocks_MP();
    benchmark_blocks_MP(p, u_overline);

    /* Benchmarking the progress made during this internship and work completed
    last year, using the same environment: 1 thread only, same A, B and p. */
    clean_file_float_integer();
    double P[5];
    P[0] = pow(2, 26) - 5;
    P[1] = pow(2, 24) - 3;
    P[2] = pow(2, 22) - 3;
    P[3] = pow(2, 20) - 3;
    P[4] = pow(2, 18) - 5;
    for (int i=0; i<5; i++){
        benchmark_float_integer(P[i], u_overline, u_b, s, t);;
    }


    return 0;
}
