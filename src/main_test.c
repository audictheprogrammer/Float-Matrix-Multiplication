#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "matrix_mul.h"
#include <assert.h>
#include <cblas.h>

// int input(int argc, char** argv){
//     if (argc != 2){
//         printf("Please input size\n");
//         char buffer[16];
//         int size;
//         if (fgets(buffer, 16, stdin) != NULL){
//             if (sscanf(buffer, "%d", &size) == 1){
//                 return size;
//             }
//         }
//     }
//
//     return atoi(argv[1]);
// }

int main(int argc, char** argv){
    const int TEST1 = 1;
    const int TEST2 = 1;
    const int TEST3 = 1;
    const int TEST4 = 1;
    const int TEST5 = 1;
    const int TEST6 = 1;
    const int TEST7 = 0;
    const int TEST8 = 0;
    const int TEST9 = 0;
    const int TEST10 = 0;


    if (TEST1){
        // Testing modular functions
        double p = pow(2, 26) - 5;

        // Precomputed constants for Modular functions
        fesetround(FE_TONEAREST);
        double u = 1.0 / p;  // Constant for SIMD
        fesetround(FE_UPWARD);
        double u_overline = 1.0 / p;  // Constant for SIMD2 and SIMD3
        fesetround(FE_TONEAREST);
        u_int32_t u_b = (int) (pow(2, 56) / p);  // Constant for Barrett

        double a = pow(2, 26);

        double SIMD1 = modulo_SIMD1(a, p, u);
        double SIMD2 = modulo_SIMD2(a, p, u_overline);
        double SIMD3 = modulo_SIMD3(a, p, u_overline);
        double Barrett = modulo_Barrett(a, p, u_b);

        printf("a = %f \n", a);
        printf("p = %f \n", p);
        printf("Correct answer is 5.0 \n");
        printf("SIMD1 returns:%f \n", SIMD1);
        printf("SIMD2 returns:%f \n", SIMD2);
        printf("SIMD3 returns:%f \n", SIMD3);
        printf("Barrett returns:%f \n\n", Barrett);

    }

    if (TEST2){
        // Testing modular functions
        double p = pow(2, 26) - 5;

        // Precomputed constants for Modular functions
        fesetround(FE_TONEAREST);
        double u = 1.0 / p;  // Constant for SIMD
        fesetround(FE_UPWARD);
        double u_overline = 1.0 / p;  // Constant for SIMD2 and SIMD3
        fesetround(FE_TONEAREST);
        u_int32_t u_b = (int) (pow(2, 56) / p);  // Constant for Barrett

        double a = (p-1) * (p-1) * 32;

        double SIMD1 = modulo_SIMD1(a, p, u);
        double SIMD2 = modulo_SIMD2(a, p, u_overline);
        double SIMD3 = modulo_SIMD3(a, p, u_overline);
        double Barrett = modulo_Barrett(a, p, u_b);

        printf("a = (p-1)(p-1) + 2^5 \n");
        printf("p = 2^26 - 5 \n");
        printf("Correct answer is 32 \n");
        printf("SIMD1 returns:%f \n", SIMD1);
        printf("SIMD2 returns:%f \n", SIMD2);
        printf("SIMD3 returns:%f \n", SIMD3);
        printf("Barrett returns:%f \n\n", Barrett);

    }

    if (TEST3){
        // Testing modular functions
        double p = pow(2, 26) - 5;

        // Precomputed constants for Modular functions
        fesetround(FE_TONEAREST);
        double u = 1.0 / p;  // Constant for SIMD
        fesetround(FE_UPWARD);
        double u_overline = 1.0 / p;  // Constant for SIMD2 and SIMD3
        fesetround(FE_TONEAREST);
        u_int32_t u_b = (int) (pow(2, 56) / p);  // Constant for Barrett

        double a = (p-1) * (p-1) * 16;  // a > 2^(25+25+4) = 2^54

        double SIMD1 = modulo_SIMD1(a, p, u);
        double SIMD2 = modulo_SIMD2(a, p, u);
        double SIMD3 = modulo_SIMD3(a, p, u);
        double Barrett = modulo_Barrett(a, p, u_b);

        printf("a = (p-1)(p-1)* 2^4 \n");
        printf("p = 2^26 - 5 \n");
        printf("Correct answer is 16 \n");
        printf("SIMD1 returns:%f \n", SIMD1);
        printf("SIMD2 returns:%f \n", SIMD2);
        printf("SIMD3 returns:%f \n", SIMD3);
        printf("Barrett returns:%f \n\n", Barrett);

    }

    if (TEST4){
        // Test for matrix multiplication
        srand(time(NULL));
        double p = pow(2, 26) - 5;

        // Precomputed constants for Modular functions
        fesetround(FE_TONEAREST);
        double u = 1.0 / p;  // Constant for SIMD
        fesetround(FE_UPWARD);
        double u_overline = 1.0 / p;  // Constant for SIMD2 and SIMD3
        u_int32_t u_b = (int) (pow(2, 56) / p);  // Constant for Barrett
        fesetround(FE_TONEAREST);

        int n = 1;

        for (int i=0; i<1; i++){

            double**A = random_matrix_2D(n, p);
            double**B = random_matrix_2D(n, p);
            double**C = zero_matrix_2D(n);  // Naive
            double**D = zero_matrix_2D(n);  // SIMD1
            double**E = zero_matrix_2D(n);  // SIMD2
            double**F = zero_matrix_2D(n);  // SIMD3
            double**G = zero_matrix_2D(n);  // Barrett


            mp_naive(A, B, C, n, p);
            mp_SIMD1(A, B, D, n, p, u);
            mp_SIMD2(A, B, E, n, p, u_overline);
            mp_SIMD3(A, B, F, n, p, u_overline);
            mp_Barrett(A, B, G, n, p, u_b);


            write_matrix(A, n, "data/Matrix_A.txt");
            write_matrix(B, n, "data/Matrix_B.txt");
            write_matrix(C, n, "data/Matrix_C.txt");  // Naive
            write_matrix(D, n, "data/Matrix_D.txt");  // SIMD1
            write_matrix(E, n, "data/Matrix_E.txt");  // SIMD2
            write_matrix(F, n, "data/Matrix_F.txt");  // SIMD3
            write_matrix(G, n, "data/Matrix_G.txt");  // Barrett

            int nb1 = equals_matrix_2D_2D(C, D, n);
            int nb2 = equals_matrix_2D_2D(C, E, n);
            int nb3 = equals_matrix_2D_2D(C, F, n);
            int nb4 = equals_matrix_2D_2D(C, G, n);

            delete_matrix_2D(&A, n);
            delete_matrix_2D(&B, n);
            delete_matrix_2D(&C, n);
            delete_matrix_2D(&D, n);
            delete_matrix_2D(&E, n);
            delete_matrix_2D(&F, n);
            delete_matrix_2D(&G, n);


            printf("i=%d \n", i);
            assert(nb1==1);
            assert(nb2==1);
            assert(nb3==1);
            assert(nb4==1);
        }

        printf("Tests passed \n");

    }

    if (TEST5){
        /* Testing the equivalence between
        if (d < 0) return d+p; else return d;
        And
        return d + (-(d<0) & (int)p)
        And
        return d + p * (d<0)
        */

        double d1 = -2;
        double d2 = 2;
        double p = 100.0;

        double res1 = d1 + ((-(d1<0)) & (u_int64_t)p);
        double res2 = d2 + ((-(d2<0)) & (u_int64_t)p);
        double res3 = d1 + p * (d1<0);
        double res4 = d2 + p * (d2<0);

        printf("res1 = %f \n", res1);
        printf("res2 = %f \n", res2);
        printf("res3 = %f \n", res3);
        printf("res4 = %f \n", res4);
    }

    if (TEST6){
        // Testing OpenBLAS's mp and my mp.
        srand(time(NULL));
        double p = pow(2, 26) - 5;
        int bitsize_p = 26;
        int n = 2;

        int blocksize = get_blocksize(bitsize_p, n);
        printf("blocksize = %d \n", blocksize);
        double** A = random_matrix_2D(n, p);
        double** B = random_matrix_2D(n, p);
        double* A_1D = convert_2D_to_1D(A, n);
        double* B_1D = convert_2D_to_1D(B, n);
        double* C = zero_matrix_1D(n*n);
        double* D = zero_matrix_1D(n*n);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,n,n,n, 1, A_1D, n, B_1D,n, 1, C,n);
        mp_ikj(A_1D, B_1D, D, n);

        print_matrix_2D(A, n);
        printf("\n");
        print_matrix_2D(B, n);
        printf("\n");
        print_matrix_1D(C, n);
        printf("\n");
        print_matrix_1D(D, n);

        printf("C == D: %d \n", equals_matrix_1D_1D(D, C, n));

        delete_matrix_2D(&A, n);
        delete_matrix_2D(&B, n);

        delete_matrix_1D(&A_1D, n);
        delete_matrix_1D(&B_1D, n);
        delete_matrix_1D(&C, n);
        delete_matrix_1D(&D, n);

    }

    if (TEST7){
        // Testing cblas_dgemm's lda ldb and ldc parameters.
        int n;

        double** A = read_matrix("data/Matrix_A_1D_2.txt", &n);
        double** B = read_matrix("data/Matrix_B_1D_2.txt", &n);
        double* A_1D = convert_2D_to_1D(A, n);
        double* B_1D = convert_2D_to_1D(B, n);
        double* C = zero_matrix_1D(n*n);
        double* D = zero_matrix_1D(n*n);
        double* E = zero_matrix_1D(n*n);

        printf("A[%d] = %f \n", n*n, A_1D[n*n]);
        printf("A[%d] = %f \n", n*n+1, A_1D[n*n+1]);
        // A_1D[n*n] = -1;
        print_matrix_1D(A_1D, n);
        print_matrix_1D(B_1D, n);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,n,n,n, 1, A_1D, 1+n, B_1D, n, 1, C,n);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,n,n,n, 1, A_1D, 2*n, B_1D, n, 1, D,n);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,n,n,n, 1, A_1D, n, B_1D, n, 1, E,n);

        printf("LDA = %d \n", n+1);
        print_matrix_1D(C, n);
        printf("LDA = %d \n", 2*n);
        print_matrix_1D(D, n);
        printf("LDA = %d \n", n);
        print_matrix_1D(E, n);


        delete_matrix_2D(&A, n);
        delete_matrix_2D(&B, n);

        delete_matrix_1D(&A_1D, n);
        delete_matrix_1D(&B_1D, n);
        delete_matrix_1D(&C, n);
        delete_matrix_1D(&D, n);
        delete_matrix_1D(&E, n);

    }
    if (TEST8){
        // Testing sub_matrix product with OpenBLAS.
        int n;
        double p = pow(2, 26) - 5;
        // double p = pow(2, 25) - 39;

        fesetround(FE_UPWARD);
        double u_overline = 1.0 / p;  // Constant for SIMD2 and SIMD3
        fesetround(FE_TONEAREST);

        double** A = read_matrix("data/Matrix_A_1D_3.txt", &n);
        double** B = read_matrix("data/Matrix_B_1D_3.txt", &n);

        int bitsize_p = get_bitsize(p);
        int b = get_blocksize(bitsize_p, n);
        printf("b = %d \n", b);

        double* A_1D = convert_2D_to_1D(A, n);
        double* B_1D = convert_2D_to_1D(B, n);
        double* C = zero_matrix_1D(n*n);
        double* D = zero_matrix_1D(n*n);
        double* E = zero_matrix_1D(n*n);

        printf("Matrix A: \n");
        print_matrix_1D(A_1D, n);
        printf("Matrix B: \n");
        print_matrix_1D(B_1D, n);

        mp_block_BLAS(A_1D, B_1D, C, n, p, u_overline, b);
        mp_block(A_1D, B_1D, D, n, p, u_overline, b);
        mp_kij(A_1D, B_1D, E, n);

        printf("Matrix C: \n");
        print_matrix_1D(C, n);
        printf("Matrix D: \n");
        print_matrix_1D(D, n);
        printf("Matrix E: \n");
        print_matrix_1D(E, n);

        delete_matrix_2D(&A, n);
        delete_matrix_2D(&B, n);

        delete_matrix_1D(&A_1D, n);
        delete_matrix_1D(&B_1D, n);
        delete_matrix_1D(&C, n);
        delete_matrix_1D(&D, n);
        delete_matrix_1D(&E, n);

    }

    if (TEST9){
        // Test for matrix product with blocks.
        srand(time(NULL));
        double p = pow(2, 26) - 5;
        // double p = pow(2, 25) - 39;


        // Precomputed constants for Modular and Blocking functions.
        fesetround(FE_UPWARD);
        double u_overline = 1.0 / p;
        int n = 256;
        int bitsize_p = get_bitsize(p);
        printf("Bitsize = %d \n", bitsize_p);
        int b = get_blocksize(bitsize_p, n);
        printf("BLocksize = %d \n", b);


        for (int i=0; i<10; i++){

            double**A = random_matrix_2D(n, p);
            double**B = random_matrix_2D(n, p);
            double* A_1D = convert_2D_to_1D(A, n);
            double* B_1D = convert_2D_to_1D(B, n);
            double*C = zero_matrix_1D(n*n);  // With BLAS
            double*D = zero_matrix_1D(n*n);  // Without BLAS
            double**E = zero_matrix_2D(n);  // Naive mp

            mp_block_BLAS(A_1D, B_1D, C, n, p, u_overline, b);
            mp_block(A_1D, B_1D, D, n, p, u_overline, b);
            mp_SIMD2(A, B, E, n, p, u_overline);


            // write_matrix(A, n, "data/Matrix_A.txt");
            // write_matrix(B, n, "data/Matrix_B.txt");
            // write_matrix(C, n, "data/Matrix_C.txt");  // Naive
            // write_matrix(D, n, "data/Matrix_D.txt");  // SIMD1
            // write_matrix(E, n, "data/Matrix_E.txt");  // SIMD2
            // write_matrix(F, n, "data/Matrix_F.txt");  // SIMD3
            // write_matrix(G, n, "data/Matrix_G.txt");  // Barrett

            int nb1 = equals_matrix_2D_1D(E, C, n);
            int nb2 = equals_matrix_2D_1D(E, D, n);

            delete_matrix_2D(&A, n);
            delete_matrix_2D(&B, n);
            delete_matrix_1D(&A_1D, n);
            delete_matrix_1D(&B_1D, n);
            delete_matrix_1D(&C, n);
            delete_matrix_1D(&D, n);
            delete_matrix_2D(&E, n);


            printf("i=%d \n", i);
            assert(nb1==1);
            assert(nb2==1);
        }

        printf("Tests passed \n");
    }

    if (TEST10){
        // Test for OpenMP
        srand(time(NULL));
        double p = pow(2, 26) - 5;

        // Precomputed constants for Modular functions
        fesetround(FE_TONEAREST);
        double u = 1.0 / p;  // Constant for SIMD
        fesetround(FE_UPWARD);
        double u_overline = 1.0 / p;  // Constant for SIMD2 and SIMD3
        u_int32_t u_b = (int) (pow(2, 56) / p);  // Constant for Barrett
        fesetround(FE_TONEAREST);

        int n = 256;

        for (int i=0; i<10; i++){

            double**A = random_matrix_2D(n, p);
            double**B = random_matrix_2D(n, p);
            double**C = zero_matrix_2D(n);  // Naive
            double**D = zero_matrix_2D(n);  // SIMD1
            double**E = zero_matrix_2D(n);  // SIMD2
            double**F = zero_matrix_2D(n);  // SIMD3
            double**G = zero_matrix_2D(n);  // Barrett


            mp_naive_MP(A, B, C, n, p);
            mp_SIMD1(A, B, D, n, p, u);
            mp_SIMD2_MP(A, B, E, n, p, u_overline);
            mp_SIMD3_MP(A, B, F, n, p, u_overline);
            mp_Barrett_MP(A, B, G, n, p, u_b);


            write_matrix(A, n, "data/Matrix_A.txt");
            write_matrix(B, n, "data/Matrix_B.txt");
            write_matrix(C, n, "data/Matrix_C.txt");  // Naive MP
            write_matrix(D, n, "data/Matrix_D.txt");  // SIMD1 MP
            write_matrix(E, n, "data/Matrix_E.txt");  // SIMD2 MP
            write_matrix(F, n, "data/Matrix_F.txt");  // SIMD3 MP
            write_matrix(G, n, "data/Matrix_G.txt");  // Barrett MP

            int nb1 = equals_matrix_2D_2D(C, D, n);
            int nb2 = equals_matrix_2D_2D(C, E, n);
            int nb3 = equals_matrix_2D_2D(C, F, n);
            int nb4 = equals_matrix_2D_2D(C, G, n);

            delete_matrix_2D(&A, n);
            delete_matrix_2D(&B, n);
            delete_matrix_2D(&C, n);
            delete_matrix_2D(&D, n);
            delete_matrix_2D(&E, n);
            delete_matrix_2D(&F, n);
            delete_matrix_2D(&G, n);


            printf("i=%d \n", i);
            assert(nb1==1);
            assert(nb2==1);
            assert(nb3==1);
            assert(nb4==1);
        }

        printf("Tests passed \n");
    }


    return 0;
}
