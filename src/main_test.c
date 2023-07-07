#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "matrix_mul.h"

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
    const int TEST4 = 1;
    // const int TEST5 = 0;
    // const int TEST6 = 0;
    // const int TEST7 = 1;


    if (TEST1){
        // Test for modulo_SIMD1
        double p = pow(2, 26) - 5;
        double u = 1.0 / p;
        double a = pow(2, 26);
        u_int32_t u_b = (int) (pow(2, 56) / p);  // Needed for Barrett

        double SIMD1 = modulo_SIMD1(a, p, u);
        double SIMD2 = modulo_SIMD2(a, p, u);
        double SIMD3 = modulo_SIMD3(a, p, u);
        double Barrett = modulo_barrett(a, p, u_b);


        printf("a = %f \n", a);
        printf("p = %f \n", p);
        printf("Correct answer is 5.0 \n");
        printf("SIMD1 returns:%f \n", SIMD1);
        printf("SIMD2 returns:%f \n", SIMD2);
        printf("SIMD3 returns:%f \n", SIMD3);
        printf("Barrett returns:%f \n\n", Barrett);


    }

    if (TEST2){
        // Test for modulo_SIMD1
        double p = pow(2, 26) - 5;
        double u = 1.0 / p;
        double a = (p-1) * (p-1) * 32;
        u_int32_t u_b = (int) (pow(2, 56) / p);

        double SIMD1 = modulo_SIMD1(a, p, u);
        double SIMD2 = modulo_SIMD2(a, p, u);
        double SIMD3 = modulo_SIMD3(a, p, u);
        double Barrett = modulo_barrett(a, p, u_b);


        printf("a = (p-1)(p-1) + 2^5 \n");
        printf("p = 2^26 - 5 \n");
        printf("Correct answer is 32 \n");
        printf("SIMD1 returns:%f \n", SIMD1);
        printf("SIMD2 returns:%f \n", SIMD2);
        printf("SIMD3 returns:%f \n", SIMD3);
        printf("Barrett returns:%f \n\n", Barrett);

    }

    if (TEST3){
        // Test for modulo_SIMD1
        double p = pow(2, 26) - 5;
        double u = 1.0 / p;
        // double a = (p-1) * (p-1) * 64;
        double a = (p-4567890);
        u_int32_t u_b = (int) (pow(2, 56) / p);

        double SIMD1 = modulo_SIMD1(a, p, u);
        double SIMD2 = modulo_SIMD2(a, p, u);
        double SIMD3 = modulo_SIMD3(a, p, u);
        double Barrett = modulo_barrett(a, p, u_b);

        printf("a = (p-1)(p-1) + 2^6 \n");
        printf("p = 2^26 - 5 \n");
        printf("Correct answer is 64 \n");
        printf("SIMD1 returns:%f \n", SIMD1);
        printf("SIMD2 returns:%f \n", SIMD2);
        printf("SIMD3 returns:%f \n", SIMD3);
        printf("Barrett returns:%f \n\n", Barrett);

    }

    if (TEST4){
        // Test for matrix multiplication
        srand(time(NULL));
        fesetround(FE_DOWNWARD);
        double p = pow(2, 26) - 5;
        double u = 1.0 / p;
        double overline_u = rint(1.0 / p);
        int n = 2;
        double**A = random_matrix(n, p);
        double**B = random_matrix(n, p);
        double**C = zero_matrix(n);  // Naive
        double**D = zero_matrix(n);  // SIMD1
        double**E = zero_matrix(n);  // SIMD2
        double**F = zero_matrix(n);  // SIMD3

        mp_naive(A, B, C, n, p);
        mp_SIMD1(A, B, D, n, p, u);
        mp_SIMD2(A, B, E, n, p, overline_u);
        mp_SIMD3(A, B, F, n, p, overline_u);

        write_matrix(A, n, "data/Matrix_A.txt");
        write_matrix(B, n, "data/Matrix_B.txt");
        write_matrix(C, n, "data/Matrix_C.txt");

        printf("Matrix A = \n");
        print_matrix(A, n);
        printf("Matrix B = \n");
        print_matrix(B, n);
        printf("Naive A*B = \n");
        print_matrix(C, n);
        printf("SIMD1 A*B = \n");
        print_matrix(D, n);
        printf("SIMD2 A*B = \n");
        print_matrix(E, n);
        printf("SIMD3 A*B = \n");
        print_matrix(D, n);
    }







    return 0;
}
