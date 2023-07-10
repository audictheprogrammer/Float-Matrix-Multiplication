#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "matrix_mul.h"
#include <assert.h>

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
    const int TEST6 = 1;
    // const int TEST7 = 1;


    if (TEST1){
        // Testing modular functions
        double p = pow(2, 26) - 5;

        // Precomputed constants for Modular functions
        fesetround(FE_DOWNWARD);
        double u = 1.0 / p;  // Constant for SIMD
        fesetround(FE_UPWARD);
        double u_overline = 1.0 / p;  // Constant for SIMD2 and SIMD3
        fesetround(FE_DOWNWARD);
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
        fesetround(FE_DOWNWARD);
        double u = 1.0 / p;  // Constant for SIMD
        fesetround(FE_UPWARD);
        double u_overline = 1.0 / p;  // Constant for SIMD2 and SIMD3
        fesetround(FE_DOWNWARD);
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
        fesetround(FE_DOWNWARD);
        double u = 1.0 / p;  // Constant for SIMD
        fesetround(FE_UPWARD);
        double u_overline = 1.0 / p;  // Constant for SIMD2 and SIMD3
        fesetround(FE_DOWNWARD);
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
        fesetround(FE_DOWNWARD);
        double u = 1.0 / p;  // Constant for SIMD
        fesetround(FE_UPWARD);
        double u_overline = 1.0 / p;  // Constant for SIMD2 and SIMD3
        fesetround(FE_DOWNWARD);
        u_int32_t u_b = (int) (pow(2, 56) / p);  // Constant for Barrett

        int n = 200;

        for (int i=0; i<10; i++){

            double**A = random_matrix(n, p);
            double**B = random_matrix(n, p);
            double**C = zero_matrix(n);  // Naive
            double**D = zero_matrix(n);  // SIMD1
            double**E = zero_matrix(n);  // SIMD2
            double**F = zero_matrix(n);  // SIMD3
            double**G = zero_matrix(n);  // Barrett


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

            int nb1 = equals_matrix(C, D, n);
            int nb2 = equals_matrix(C, E, n);
            int nb3 = equals_matrix(C, F, n);
            int nb4 = equals_matrix(C, G, n);

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
        return d + p()
        */

        double d1 = -2;
        double d2 = 2;
        double p = 100.0;

        double a = 1024.5;
        double b = 0b1111111111111111111111111111111111111111111111111111111111111111;
        double temp = b && a;

        // 0b 1111 1111 ... 1111 for double = ...

        double res1 = d1 + ((-(d1<0)) && p);
        double res2 = d2 + ((-(d2<0)) && p);

        printf("res1 = %f \n", res1);
        printf("res2 = %f \n", res2);
        printf("b = %f \n", b);
        printf("temp = %f \n", temp);
    }

    if (TEST6){
        /* Testing FE_DOWNWARD and FE_UPWARD */
        double p = pow(2, 26) - 5;
        // double u1 = (p-1) * (p-1) * (p-1);
        double u1 = -(p-1) * (p-1) * (p-1) * 64;
        fesetround(FE_DOWNWARD);
        double u2 = -(p-1) * (p-1) * (p-1) * 64;
        fesetround(FE_UPWARD);
        double u3 = -(p-1) * (p-1) * (p-1) * 64;
        fesetround(FE_TOWARDZERO);
        double u4 = -(p-1) * (p-1) * (p-1) * 64;
        fesetround(FE_TONEAREST);
        double u5 = -(p-1) * (p-1) * (p-1) * 64;

        printf("u1=%f \n", u1);
        printf("u2=%f \n", u2);
        printf("u3=%f \n", u3);
        printf("u4=%f \n", u3);
        printf("u5=%f \n", u3);


    }

    return 0;
}
