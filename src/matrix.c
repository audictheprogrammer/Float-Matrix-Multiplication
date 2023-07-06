#include "matrix.h"
#include "matrix_mul.h"

// Basic operations
double** zero_matrix(int size){
    /* Create a square matrix of size "size" */
    double** mat = (double**) malloc(sizeof(double*) * size);
    if (mat == NULL){
        printf("Malloc error !\n");
        return NULL;
    }
    for (int i = 0; i < size; i++){
        mat[i] = (double*) malloc(sizeof(double) * size);
        if (mat[i] == NULL){
            printf("Malloc error !\n");
            free(mat);
            return NULL;
        }
        for (int j = 0; j < size; j++){
            mat[i][j] = 0;
        }
    }

    return mat;
}


double** random_matrix(int size, double p){
    /* Create a square matrix of size "size" and fill
    with numbers between 0 and p-1 */

    // Create the matrix
    double** mat = zero_matrix(size);
    // Fill the matrix
    if (mat == NULL){
        return NULL;
    }
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            // int rand_int = rand();
            // mat[i][j] = modulo_SIMD1(rand_int, p);
            mat[i][j] = (double) (rand() % ((int) p));

            // mat[i][j] = (double)rand()/(double)(RAND_MAX/p);


            // mat[i][j] = rand() % mod;
        }
    }

    return mat;
}

void delete_matrix(double*** mat, int size){
    /* Delete the square matrix mat and set mat to NULL*/
    if (*mat == NULL){
        printf("Matrix is NULL cannot delete !\n");
        return;
    }
    for (int i = 0; i < size; i++){
        free((*mat)[i]);
    }
    free(*mat);
    *mat = NULL;
}


double** transpose_matrix(double** mat, int size){
    /* Tranpose the square matrix mat*/
    double** mat_t = zero_matrix(size);
    if (mat_t == NULL){
        printf("Cannot transpose the matrix !\n");
    }
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            mat_t[i][j] = mat[j][i];
        }
    }

    return mat_t;

}

void print_matrix(double** mat, int size){
    /* Print the matrix. Here is an example:
    [ 1 2 3 ]
    [ 4 5 6 ]
    [ 7 8 9 ]
    */
    if (mat == NULL){
        printf("[ ]\n");
        return;
    }

    for (int i = 0; i < size; i++){
        printf("[ ");
        for (int j = 0; j < size; j++){
            printf("%f ", mat[i][j]);
        }
        printf("]\n");
    }

    return ;
}

void write_matrix(double** mat, int size, char* filename){
    /* Write the matrix into a file. Here is an example:
    2
    [ 1 2 ]
    [ 3 4 ]
    */
    FILE* f = fopen(filename, "w");
    if (f == NULL){
        printf("Filename is incorrect !\n");
        return;
    }

    fprintf(f, "%d\n", size);
    for (int i = 0; i < size; i++){
            fprintf(f, "[ ");
        for(int j = 0; j < size; j++){
            fprintf(f, "%f ", mat[i][j]);
        }
        fprintf(f, "]\n");
    }

    fclose(f);
    return;
}

double** read_matrix(char* filename, int* size){
    /* Read a file, modify size value and return the matrix */
    char buffer[1024];
    FILE* f = fopen(filename, "r");
    if (f == NULL){
        return NULL;
    }

    // Treating first line.
    if (fgets(buffer, 1025, f) == NULL){
        printf("fgets error !\n");
        return NULL;
    }
    if (sscanf(buffer, "%d", size) != 1){
        printf("sscanf error !\n");
        return NULL;
    }
    double** mat = zero_matrix(*size);
    if (mat == NULL){
        return NULL;
    }

    // Treating next lines.
    for (int i = 0; i < *size; i++){
        fgets(buffer, 1025, f);
        char* strToken = strtok(buffer, "[] \n");
        int j = 0;
        while (strToken) {
            sscanf(strToken, "%lf", &(mat[i][j]));
            strToken = strtok(NULL, " ");
            j++;
        }
    }

    return mat;
}

int equals_matrix(double** A, double** B, int size){

    /* Check if A equals B. The return values are:
    0: different
    1: equals
    */
    if (A == NULL && B == NULL){
        return 1;
    }
    if (A == NULL || B == NULL){
        return 0;
    }

    for (int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            if (A[i][j] != B[i][j]){
                printf("A[%d][%d] = %f \n", i, j, A[i][j]);
                printf("B[%d][%d] = %f \n", i, j, B[i][j]);
                return 0;
            }
        }
    }

    return 1;
}

int equals_matrix_file(char* filename1, char* filename2){
    /* Check if A equals B. The return values are:
    0: different
    1: equals
    */
    int size1;
    int size2;
    double** mat1 = read_matrix(filename1, &size1);
    double** mat2 = read_matrix(filename2, &size2);

    return equals_matrix(mat1, mat2, size1);
}
