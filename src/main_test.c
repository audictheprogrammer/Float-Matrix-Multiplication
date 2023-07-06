#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "matrix_mul.h"

int input(int argc, char** argv){
    if (argc != 2){
        printf("Please input size\n");
        char buffer[16];
        int size;
        if (fgets(buffer, 16, stdin) != NULL){
            if (sscanf(buffer, "%d", &size) == 1){
                return size;
            }
        }
    }

    return atoi(argv[1]);
}

int main(int argc, char** argv){
    const int TEST1 = 0;
    // const int TEST2 = 0;
    // const int TEST3 = 0;
    // const int TEST4 = 0;
    // const int TEST5 = 0;
    // const int TEST6 = 0;
    // const int TEST7 = 1;


    if (TEST1) {
        srand(time(NULL));
        u_int32_t p = 1073741827;
        int size = input(argc, argv);

        u_int64_t** mat = create_matrix_mod(size, p);
        print_matrix(mat, size);

        u_int64_t** mat_t = transpose_matrix(mat, size);
        print_matrix(mat_t, size);

        delete_matrix(&mat, size);
        delete_matrix(&mat_t, size);
    }


    return 0;
}
