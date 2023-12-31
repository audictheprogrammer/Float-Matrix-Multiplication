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
    const int TEST1 = 0;
    const int TEST2 = 0;
    const int TEST3 = 0;
    const int TEST4 = 0;
    const int TEST5 = 0;
    const int TEST6 = 0;
    const int TEST7 = 0;
    const int TEST8 = 0;
    const int TEST9 = 1;
    const int TEST10 = 1;
    const int TEST11 = 1;


    if (TEST1){
        // Testing modular functions
        // Barrett works with these 5 p. Test passed.

        // double p = pow(2, 26) - 5;
        // double p = pow(2, 24) - 3;
        // double p = pow(2, 22) - 3;
        double p = pow(2, 20) - 3;
        // double p = pow(2, 18) - 5;


        // Precomputed constants for Modular functions
        double u = 1.0 / p;  // Constant for SIMD
        fesetround(FE_UPWARD);
        double u_overline = 1.0 / p;  // Constant for SIMD2 and SIMD3
        fesetround(FE_TONEAREST);
        // Constants for Barrett which is the integer version of SIMD
        u_int32_t t = 32;
        u_int32_t s = 30 + get_bitsize(p) - t;
        u_int32_t u_b = (int) (pow(2, s+t) / p);

        printf("s = %d, t = %d \n", s, t);
        double A[5];
        A[0] = 150007655597277077;  // 2pow58 bitsize(p) >= 26
        A[1] = 58528953432971872;  // 2pow56 bitsize(p) >= 24
        A[2] = 10490987492010470;  // 2pow54 bitsize(p) >= 22
        A[3] = 3088868763674166;  // 2pow52 bitsize(p) >= 20
        A[4] = 613220816403662;  // 2pow50 bitsize(p) >= 18

        for (int i=0; i<5; i++){
            double SIMD1 = modulo_SIMD1(A[i], p, u);
            double SIMD2 = modulo_SIMD2(A[i], p, u_overline);
            double SIMD3 = modulo_SIMD3(A[i], p, u_overline);
            double Barrett = modulo_Barrett(A[i], p, u_b, s, t);
            printf("%f%%%f \n", A[i], p);
            printf("SIMD1 returns: %f \n", SIMD1);
            printf("SIMD2 returns: %f \n", SIMD2);
            printf("SIMD3 returns: %f \n", SIMD3);
            printf("Barrett returns: %f \n\n", Barrett);
        }

    }

    if (TEST2){
        // Testing modular functions
        double p = pow(2, 26) - 5;

        // Precomputed constants for Modular functions
        double u = 1.0 / p;  // Constant for SIMD
        fesetround(FE_UPWARD);
        double u_overline = 1.0 / p;  // Constant for SIMD2 and SIMD3
        fesetround(FE_TONEAREST);
        // Constants for Barrett which is the integer version of SIMD
        u_int32_t t = 32;
        u_int32_t s = 30 + get_bitsize(p) - t;
        u_int32_t u_b = (int) (pow(2, s+t) / p);

        double a = (p-1) * (p-1) * 32;

        double SIMD1 = modulo_SIMD1(a, p, u);
        double SIMD2 = modulo_SIMD2(a, p, u_overline);
        double SIMD3 = modulo_SIMD3(a, p, u_overline);
        double Barrett = modulo_Barrett(a, p, u_b, s, t);

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
        double u = 1.0 / p;  // Constant for SIMD
        fesetround(FE_UPWARD);
        double u_overline = 1.0 / p;  // Constant for SIMD2 and SIMD3
        fesetround(FE_TONEAREST);
        // Constants for Barrett which is the integer version of SIMD
        u_int32_t t = 32;
        u_int32_t s = 30 + get_bitsize(p) - t;
        u_int32_t u_b = (int) (pow(2, s+t) / p);

        double a = (p-1) * (p-1) * 16;  // a > 2^(25+25+4) = 2^54

        double SIMD1 = modulo_SIMD1(a, p, u);
        double SIMD2 = modulo_SIMD2(a, p, u);
        double SIMD3 = modulo_SIMD3(a, p, u);
        double Barrett = modulo_Barrett(a, p, u_b, s, t);

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
        double u = 1.0 / p;  // Constant for SIMD
        fesetround(FE_UPWARD);
        double u_overline = 1.0 / p;  // Constant for SIMD2 and SIMD3
        fesetround(FE_TONEAREST);
        // Constants for Barrett which is the integer version of SIMD
        u_int32_t t = 32;
        u_int32_t s = 30 + get_bitsize(p) - t;
        u_int32_t u_b = (int) (pow(2, s+t) / p);

        int n = 1;

        for (int i=0; i<1; i++){

            double* A = random_matrix_1D(n, p);
            double* B = random_matrix_1D(n, p);
            double* C = zero_matrix_1D(n);  // Naive
            double* D = zero_matrix_1D(n);  // SIMD1
            double* E = zero_matrix_1D(n);  // SIMD2
            double* F = zero_matrix_1D(n);  // SIMD3
            double* G = zero_matrix_1D(n);  // Barrett


            mp_naive(A, B, C, n, p);
            mp_SIMD1(A, B, D, n, p, u);
            mp_SIMD2(A, B, E, n, p, u_overline);
            mp_SIMD3(A, B, F, n, p, u_overline);
            mp_Barrett(A, B, G, n, p, u_b, s, t);


            write_matrix_1D(A, n, "data/Matrix_A.txt");
            write_matrix_1D(B, n, "data/Matrix_B.txt");
            write_matrix_1D(C, n, "data/Matrix_C.txt");  // Naive
            write_matrix_1D(D, n, "data/Matrix_D.txt");  // SIMD1
            write_matrix_1D(E, n, "data/Matrix_E.txt");  // SIMD2
            write_matrix_1D(F, n, "data/Matrix_F.txt");  // SIMD3
            write_matrix_1D(G, n, "data/Matrix_G.txt");  // Barrett

            int nb1 = equals_matrix_1D_1D(C, D, n);
            int nb2 = equals_matrix_1D_1D(C, E, n);
            int nb3 = equals_matrix_1D_1D(C, F, n);
            int nb4 = equals_matrix_1D_1D(C, G, n);

            delete_matrix_1D(&A, n);
            delete_matrix_1D(&B, n);
            delete_matrix_1D(&C, n);
            delete_matrix_1D(&D, n);
            delete_matrix_1D(&E, n);
            delete_matrix_1D(&F, n);
            delete_matrix_1D(&G, n);


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
        // Testing OpenBLAS's mp.
        srand(time(NULL));
        double p = pow(2, 26) - 5;
        int bitsize_p = 26;
        int n = 2;

        int blocksize = get_blocksize(bitsize_p, n);
        printf("blocksize = %d \n", blocksize);
        double* A = random_matrix_1D(n, p);
        double* B = random_matrix_1D(n, p);
        double* C = zero_matrix_1D(n);
        double* D = zero_matrix_1D(n);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,n,n,n, 1, A, n, B,n, 1, C,n);
        mp_ikj(A, B, D, n);

        print_matrix_1D(A, n);
        printf("\n");
        print_matrix_1D(B, n);
        printf("\n");
        print_matrix_1D(C, n);
        printf("\n");
        print_matrix_1D(D, n);

        printf("C == D: %d \n", equals_matrix_1D_1D(D, C, n));

        delete_matrix_1D(&A, n);
        delete_matrix_1D(&B, n);
        delete_matrix_1D(&C, n);
        delete_matrix_1D(&D, n);

    }

    if (TEST7){
        // Testing cblas_dgemm's lda ldb and ldc parameters.
        int n;
        double* A = read_matrix_1D("data/Matrix_A_2.txt", &n);
        double* B = read_matrix_1D("data/Matrix_B_2.txt", &n);
        if (A == NULL || B == NULL)
            return 1;
        double* C = zero_matrix_1D(n);
        double* D = zero_matrix_1D(n);
        double* E = zero_matrix_1D(n);

        // printf("A[%d] = %f \n", n*n, A[n*n]);
        // printf("A[%d] = %f \n", n*n+1, A[n*n+1]);
        // A[n*n] = -1;
        print_matrix_1D(A, n);
        print_matrix_1D(B, n);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,n,n,n, 1, A, 1+n, B, n, 1, C,n);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,n,n,n, 1, A, 2*n, B, n, 1, D,n);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,n,n,n, 1, A, n, B, n, 1, E,n);

        printf("LDA = %d \n", n+1);
        print_matrix_1D(C, n);
        printf("LDA = %d \n", 2*n);
        print_matrix_1D(D, n);
        printf("LDA = %d \n", n);
        print_matrix_1D(E, n);


        delete_matrix_1D(&A, n);
        delete_matrix_1D(&B, n);
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

        double* A = read_matrix_1D("data/Matrix_A_3.txt", &n);
        double* B = read_matrix_1D("data/Matrix_B_3.txt", &n);

        int bitsize_p = get_bitsize(p);
        int b = get_blocksize(bitsize_p, n);
        printf("b = %d \n", b);

        double* C = zero_matrix_1D(n);
        double* D = zero_matrix_1D(n);
        double* E = zero_matrix_1D(n);

        printf("Matrix A: \n");
        print_matrix_1D(A, n);
        printf("Matrix B: \n");
        print_matrix_1D(B, n);

        mp_block_BLAS(A, B, C, n, p, u_overline, b);
        mp_block(A, B, D, n, p, u_overline, b);
        mp_kij(A, B, E, n);

        printf("Matrix C: \n");
        print_matrix_1D(C, n);
        printf("Matrix D: \n");
        print_matrix_1D(D, n);
        printf("Matrix E: \n");
        print_matrix_1D(E, n);

        delete_matrix_1D(&A, n);
        delete_matrix_1D(&B, n);
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
        fesetround(FE_TONEAREST);
        int n = 256;
        int bitsize_p = get_bitsize(p);
        printf("Bitsize = %d \n", bitsize_p);
        int b = get_blocksize(bitsize_p, n);
        printf("BLocksize = %d \n", b);


        for (int i=0; i<10; i++){

            double* A = random_matrix_1D(n, p);
            double* B = random_matrix_1D(n, p);
            double* C = zero_matrix_1D(n);  // With BLAS
            double* D = zero_matrix_1D(n);  // Without BLAS
            double* E = zero_matrix_1D(n);  // Naive mp

            mp_block_BLAS(A, B, C, n, p, u_overline, b);
            mp_block(A, B, D, n, p, u_overline, b);
            mp_SIMD2(A, B, E, n, p, u_overline);


            // write_matrix_1D(A, n, "data/Matrix_A.txt");
            // write_matrix_1D(B, n, "data/Matrix_B.txt");
            // write_matrix_1D(C, n, "data/Matrix_C.txt");  // Naive
            // write_matrix_1D(D, n, "data/Matrix_D.txt");  // SIMD1
            // write_matrix_1D(E, n, "data/Matrix_E.txt");  // SIMD2
            // write_matrix_1D(F, n, "data/Matrix_F.txt");  // SIMD3
            // write_matrix_1D(G, n, "data/Matrix_G.txt");  // Barrett

            int nb1 = equals_matrix_1D_1D(E, C, n);
            int nb2 = equals_matrix_1D_1D(E, D, n);

            delete_matrix_1D(&A, n);
            delete_matrix_1D(&B, n);
            delete_matrix_1D(&C, n);
            delete_matrix_1D(&D, n);
            delete_matrix_1D(&E, n);


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
        double u = 1.0 / p;  // Constant for SIMD
        fesetround(FE_UPWARD);
        double u_overline = 1.0 / p;  // Constant for SIMD2 and SIMD3
        fesetround(FE_TONEAREST);
        // Constants for Barrett which is the integer version of SIMD
        u_int32_t t = 32;
        u_int32_t s = 30 + get_bitsize(p) - t;
        u_int32_t u_b = (int) (pow(2, s+t) / p);

        int n = 128;

        for (int i=0; i<10; i++){

            double* A = random_matrix_1D(n, p);
            double* B = random_matrix_1D(n, p);
            double* C = zero_matrix_1D(n);  // Naive
            double* D = zero_matrix_1D(n);  // SIMD1
            double* E = zero_matrix_1D(n);  // SIMD2
            double* F = zero_matrix_1D(n);  // SIMD3
            double* G = zero_matrix_1D(n);  // Barrett


            mp_naive_MP(A, B, C, n, p);
            mp_SIMD1(A, B, D, n, p, u);
            mp_SIMD2_MP(A, B, E, n, p, u_overline);
            mp_SIMD3_MP(A, B, F, n, p, u_overline);
            mp_Barrett_MP(A, B, G, n, p, u_b, s, t);


            write_matrix_1D(A, n, "data/Matrix_A.txt");
            write_matrix_1D(B, n, "data/Matrix_B.txt");
            write_matrix_1D(C, n, "data/Matrix_C.txt");  // Naive MP
            write_matrix_1D(D, n, "data/Matrix_D.txt");  // SIMD1 MP
            write_matrix_1D(E, n, "data/Matrix_E.txt");  // SIMD2 MP
            write_matrix_1D(F, n, "data/Matrix_F.txt");  // SIMD3 MP
            write_matrix_1D(G, n, "data/Matrix_G.txt");  // Barrett MP

            int nb1 = equals_matrix_1D_1D(C, D, n);
            int nb2 = equals_matrix_1D_1D(C, E, n);
            int nb3 = equals_matrix_1D_1D(C, F, n);
            int nb4 = equals_matrix_1D_1D(C, G, n);

            delete_matrix_1D(&A, n);
            delete_matrix_1D(&B, n);
            delete_matrix_1D(&C, n);
            delete_matrix_1D(&D, n);
            delete_matrix_1D(&E, n);
            delete_matrix_1D(&F, n);
            delete_matrix_1D(&G, n);


            printf("i=%d \n", i);
            assert(nb1==1);
            assert(nb2==1);
            assert(nb3==1);
            assert(nb4==1);
        }

        printf("Tests passed \n");
    }

    if (TEST11){
        // Testing integer mp and float mp.

        srand(time(NULL));
        // double p = 94906249;  // 2^{26} < p < 2^{26.5}, please set b = 1
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
        // Constants for Barrett which is the integer version of SIMD
        u_int32_t t = 32;
        u_int32_t s = 30 + get_bitsize(p) - t;
        u_int32_t u_b = (int) (pow(2, s+t) / p);

        printf("s = %d, t = %d \n", s, t);

        openblas_set_num_threads(1);

        int n = 512;
        int bitsize_p = get_bitsize(p);
        int b = get_blocksize(bitsize_p, n);

        double* A = random_matrix_1D(n, p);
        double* B = random_matrix_1D(n, p);

        u_int64_t* A_int = convert_float_to_integer(A, n);
        u_int64_t* B_temp = convert_float_to_integer(B, n);
        u_int64_t* B_int = transpose_matrix_integer(B_temp, n);

        double* C = zero_matrix_1D(n);
        u_int64_t* D = zero_matrix_1D_integer(n);

        // PRODUCT STARTS **********************************
        mp_float(A, B, C, n, p, u_overline, b);
        mp_integer(A_int, B_int, D, n, p, u_b, s, t);
        // PRODUCT ENDS ************************************

        // printf("Matrix A: \n");
        // print_matrix_1D_integer(A_int, n);
        // printf("Matrix B: \n");
        // print_matrix_1D_integer(B_int, n);
        //
        // printf("Matrix C: \n");
        // print_matrix_1D(C, n);
        // printf("Matrix D: \n");
        // print_matrix_1D_integer(D, n);

        write_matrix_1D_integer(A_int, n, "data/Matrix_A.txt");
        write_matrix_1D_integer(B_int, n, "data/Matrix_B.txt");
        write_matrix_1D(C, n, "data/Matrix_C.txt");  // Double
        write_matrix_1D_integer(D, n, "data/Matrix_D.txt");  // Integer


        int nb1 = equals_matrix_float_integer(C, D, n);

        delete_matrix_1D(&A, n);
        delete_matrix_1D_integer(&A_int, n);
        delete_matrix_1D(&B, n);
        delete_matrix_1D_integer(&B_temp, n);
        delete_matrix_1D_integer(&B_int, n);
        delete_matrix_1D(&C, n);
        delete_matrix_1D_integer(&D, n);

        assert(nb1==1);
        printf("Test passed \n");

    }


    return 0;
}
