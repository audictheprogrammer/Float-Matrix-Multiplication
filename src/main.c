#include "matrix.h"
#include "matrix_mul.h"
#include <math.h>


// Matrix Product with different order of loops.
double benchmark_ijk(double** A, double** B, int n){
    double** C = zero_matrix_2D(n);

    clock_t initial = clock();
    mp_ijk(A, B, C, n);
    clock_t final = clock();

    delete_matrix_2D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
}

double benchmark_kij(double** A, double** B, int n){
    double** C = zero_matrix_2D(n);

    clock_t initial = clock();
    mp_kij(A, B, C, n);
    clock_t final = clock();

    delete_matrix_2D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
}

double benchmark_jki(double** A, double** B, int n){
    double** C = zero_matrix_2D(n);

    clock_t initial = clock();
    mp_jki(A, B, C, n);
    clock_t final = clock();

    delete_matrix_2D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
}

double benchmark_ikj(double** A, double** B, int n){
    double** C = zero_matrix_2D(n);

    clock_t initial = clock();
    mp_ikj(A, B, C, n);
    clock_t final = clock();

    delete_matrix_2D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
}

double benchmark_jik(double** A, double** B, int n){
    double** C = zero_matrix_2D(n);

    clock_t initial = clock();
    mp_jik(A, B, C, n);
    clock_t final = clock();

    delete_matrix_2D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
}

double benchmark_kji(double** A, double** B, int n){
    double** C = zero_matrix_2D(n);

    clock_t initial = clock();
    mp_kji(A, B, C, n);
    clock_t final = clock();

    delete_matrix_2D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
}


// Matrix Product with modulos-.

double benchmark_mod_naive(double** A, double** B, int n, double p){
    double** C = zero_matrix_2D(n);

    clock_t initial = clock();
    mp_naive(A, B, C, n, p);
    clock_t final = clock();

    delete_matrix_2D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
}

double benchmark_mod_SIMD1(double** A, double** B, int n, double p, double u){
    double** C = zero_matrix_2D(n);

    clock_t initial = clock();
    mp_SIMD1(A, B, C, n, p, u);
    clock_t final = clock();

    delete_matrix_2D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
}

double benchmark_mod_SIMD2(double** A, double** B, int n, double p, double u){
    double** C = zero_matrix_2D(n);

    clock_t initial = clock();
    mp_SIMD2(A, B, C, n, p, u);
    clock_t final = clock();

    delete_matrix_2D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
}

 double benchmark_mod_SIMD3(double** A, double** B, int n, double p, double u){
    double** C = zero_matrix_2D(n);

    clock_t initial = clock();
    mp_SIMD3(A, B, C, n, p, u);
    clock_t final = clock();

    delete_matrix_2D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
}

double benchmark_mod_Barrett(double** A, double** B, int n, double p, double u){
    double** C = zero_matrix_2D(n);

    clock_t initial = clock();
    mp_Barrett(A, B, C, n, p, u);
    clock_t final = clock();

    delete_matrix_2D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
}

// Matrix Product with Blocks.

double benchmark_blocks_NoBLAS(double* A, double* B, int n, double p, double u, int b){
    double* C = zero_matrix_1D(n*n);

    clock_t initial = clock();
    mp_block(A, B, C, n, p, u, b);
    clock_t final = clock();

    delete_matrix_1D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
}

double benchmark_blocks_BLAS(double* A, double* B, int n, double p, double u, int b){
    double* C = zero_matrix_1D(n*n);

    clock_t initial = clock();
    mp_block_BLAS(A, B, C, n, p, u, b);
    clock_t final = clock();

    delete_matrix_1D(&C, n);
    return ((double) (final - initial)) / CLOCKS_PER_SEC;
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
    int m = 5;  // Executes m times each algo
    for (int i=8; i<11; i++){
        int n = (int) pow(2, i);
        double sum_ijk = 0;
        double sum_kij = 0;
        double sum_jki = 0;
        double sum_ikj = 0;
        double sum_jik = 0;
        double sum_kji = 0;

        for (int j=0; j<m; j++){
            double**A = random_matrix_2D(n, p);
            double**B = random_matrix_2D(n, p);
            sum_ijk += benchmark_ijk(A, B, n);
            sum_kij += benchmark_kij(A, B, n);
            sum_jki += benchmark_jki(A, B, n);
            sum_ikj += benchmark_ikj(A, B, n);
            sum_jik += benchmark_jik(A, B, n);
            sum_kji += benchmark_kji(A, B, n);

            delete_matrix_2D(&A, n);
            delete_matrix_2D(&B, n);
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

void benchmark_modulos(double p, double u, double u_overline, double u_b){
    /* Benchmarking different modulos.
    The most efficient one SIMD2.
    */
    int m = 5;  // Executes m times each algo
    for (int i=8; i<11; i++){
        int n = (int) pow(2, i);
        double sum_mod_naive = 0;
        double sum_mod_SIMD1 = 0;
        double sum_mod_SIMD2 = 0;
        double sum_mod_SIMD3 = 0;
        double sum_mod_Barrett = 0;

        for (int j=0; j<m; j++){
            double**A = random_matrix_2D(n, p);
            double**B = random_matrix_2D(n, p);
            sum_mod_naive += benchmark_mod_naive(A, B, n, p);
            sum_mod_SIMD1 += benchmark_mod_SIMD1(A, B, n, p, u);
            sum_mod_SIMD2 += benchmark_mod_SIMD2(A, B, n, p, u_overline);
            sum_mod_SIMD3 += benchmark_mod_SIMD3(A, B, n, p, u_overline);
            sum_mod_Barrett += benchmark_mod_Barrett(A, B, n, p, u_b);

            delete_matrix_2D(&A, n);
            delete_matrix_2D(&B, n);
        }

        printf("\n");
        write_benchmark_time("data/benchmark_modulo_naive.txt", "Mod Naive", n, sum_mod_naive/m);
        write_benchmark_time("data/benchmark_modulo_SIMD1.txt", "Mod SIMD1", n, sum_mod_SIMD1/m);
        write_benchmark_time("data/benchmark_modulo_SIMD2.txt", "Mod SIMD2", n, sum_mod_SIMD2/m);
        write_benchmark_time("data/benchmark_modulo_SIMD3.txt", "Mod SIMD3", n, sum_mod_SIMD3/m);
        write_benchmark_time("data/benchmark_modulo_Barrett.txt", "Mod Barrett", n, sum_mod_Barrett/m);

    }
}

void benchmark_blocks(double p, double u_overline){
    /* Benchmarking different modulos.
    */
    int m = 5;  // Executes m times each algo
    for (int i=8; i<11; i++){
        int n = (int) pow(2, i);
        int b = get_blocksize(get_bitsize(p), n);

        double sum_blocks = 0;
        double sum_blocks_BLAS = 0;

        for (int j=0; j<m; j++){
            double*A = random_matrix_1D(n, p);
            double*B = random_matrix_1D(n, p);
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

void clean_file_blocks(){
    char noms[2][64] = {"data/benchmark_blocks_BLAS.txt", "data/benchmark_blocks_NoBLAS.txt"};

    for (int i=0; i<2; i++){
        FILE *f = fopen(noms[i], "w");
        fclose(f);
    }

}


int main(){
    // Initialization
    srand(time(NULL));
    double p = pow(2, 26) - 5;

    // Precomputed constants for Modular functions
    double u = 1.0 / p;  // Constant for SIMD
    fesetround(FE_UPWARD);
    double u_overline = 1.0 / p;  // Constant for SIMD2 and SIMD3
    fesetround(FE_TONEAREST);
    u_int32_t u_b = (int) (pow(2, 56) / p);  // Constant for Barrett


    // Benchmarking order of loop.
    // 07/07/23 13:27 I did a benchmark for 5
    // clean_file_loops();
    // benchmark_loops_order(p);

    // Benchmarking different modulos.
    //
    // clean_file_modulos();
    // benchmark_modulos(p, u, u_overline, u_b);

    // Benchmarking blocks.
    clean_file_blocks();
    benchmark_blocks(p, u_overline);


    return 0;
}
