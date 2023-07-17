#include "matrix.h"
#include "matrix_mul.h"

// Basic operations
double** zero_matrix_2D(int n){
    /* Allocate a 2 dimensional square matrix. */
    double** mat = (double**) malloc(sizeof(double*) * n);
    if (mat == NULL){
        printf("Malloc error !\n");
        return NULL;
    }
    for (int i = 0; i < n; i++){
        mat[i] = (double*) malloc(sizeof(double) * n);
        if (mat[i] == NULL){
            printf("Malloc error !\n");
            free(mat);
            return NULL;
        }
        for (int j = 0; j < n; j++){
            mat[i][j] = 0;
        }
    }

    return mat;
}

double* zero_matrix_1D(int n){
    /* Allocate a 1 dimensional square matrix. */
    double* mat = (double*) malloc(sizeof(double) * n);
    if (mat == NULL){
        printf("Malloc error !\n");
        return NULL;
    }
    for (int i = 0; i < n; i++){
        mat[i] = 0;
    }

    return mat;
}


double** random_matrix_2D(int n, double p){
    /* Create a 2 dimensional square matrix and fill
    with numbers between 0 and p-1. */

    double** mat = zero_matrix_2D(n);
    // ALlocate the matrix
    // Fill the matrix
    if (mat == NULL){
        return NULL;
    }
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            mat[i][j] = (double) (rand() % ((int) p));
        }
    }

    return mat;
}

double* random_matrix_1D(int n, double p){
    /* Create a 1 dimensional square matrix and fill
    with numbers between 0 and p-1. */

    // Allocate the matrix
    double* mat = zero_matrix_1D(n);
    // Fill the matrix
    if (mat == NULL){
        return NULL;
    }
    for (int i = 0; i < n*n; i++){
            mat[i] = (double) (rand() % ((int) p));
        }

    return mat;
}


double** convert_1D_to_2D(double* A, int n){
    /* Converts a 1 dimensional square matrix to
    a 2 dimensional matrix. */
    double** mat = zero_matrix_2D(n);
    int i = 0;
    int j = 0;

    for (int k = 0; k < n*n; k++){
        mat[i][j] = A[k];
        if (j+1 == n){
            j = 0;
            i++;
        } else {
            j++;
        }
    }
    return mat;
}

double* convert_2D_to_1D(double** A, int n){
    /* Converts a 2 dimensional square matrix to
    a 1 dimensional matrix. */
    double* mat = zero_matrix_1D(n*n);
    int count = 0;
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            mat[count] = A[i][j];
            count++;
        }
    }
    return mat;
}


void delete_matrix_2D(double*** mat, int n){
    /* Delete the 2D matrix mat and set it to NULL*/
    if (*mat == NULL){
        printf("Matrix is NULL cannot delete !\n");
        return;
    }
    for (int i = 0; i < n; i++){
        free((*mat)[i]);
    }
    free(*mat);
    *mat = NULL;
}

void delete_matrix_1D(double** mat, int n){
    /* Delete the 1D matrix and set it to NULL*/
    if (*mat == NULL){
        printf("Matrix is NULL cannot delete !\n");
        return;
    }
    free(*mat);
    *mat = NULL;
}


double** transpose_matrix(double** mat, int n){
    /* Tranpose the square matrix mat*/
    double** mat_t = zero_matrix_2D(n);
    if (mat_t == NULL){
        printf("Cannot transpose the matrix !\n");
    }
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            mat_t[i][j] = mat[j][i];
        }
    }
    return mat_t;
}


void print_matrix_2D(double** mat, int n){
    /* Print the matrix. Here is an example:
    [ 1 2 3 ]
    [ 4 5 6 ]
    [ 7 8 9 ]
    */
    if (mat == NULL){
        printf("[ ]\n");
        return;
    }
    for (int i = 0; i < n; i++){
        printf("[ ");
        for (int j = 0; j < n; j++){
            printf("%f ", mat[i][j]);
        }
        printf("]\n");
    }
    return ;
}


void print_matrix_1D(double* mat, int n){
    /* Print the matrix. Here is an example:
    [ 1 2 3 ]
    [ 4 5 6 ]
    [ 7 8 9 ]
    */
    if (mat == NULL){
        printf("[ ]\n");
        return;
    }
    for (int i = 0; i < n*n; i++){
        if (i % n == 0){
            printf("[ ");
        }
        printf("%f ", mat[i]);
        if ((i+1) % n == 0){
            printf("]\n");
        }
    }
    return ;
}

void write_matrix(double** mat, int n, char* filename){
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

    fprintf(f, "%d\n", n);
    for (int i = 0; i < n; i++){
            fprintf(f, "[ ");
        for(int j = 0; j < n; j++){
            fprintf(f, "%f ", mat[i][j]);
        }
        fprintf(f, "]\n");
    }

    fclose(f);
    return;
}

double** read_matrix(char* filename, int* n){
    /* Read the file, update n, then return the matrix. */
    char buffer[1024];
    FILE* f = fopen(filename, "r");
    if (f == NULL){
        return NULL;
    }

    // Treating first line.
    if (fgets(buffer, 1025, f) == NULL){
        printf("fgets error !\n");
        fclose(f);
        return NULL;
    }
    if (sscanf(buffer, "%d", n) != 1){
        printf("sscanf error !\n");
        fclose(f);
        return NULL;
    }
    double** mat = zero_matrix_2D(*n);
    if (mat == NULL){
        fclose(f);
        return NULL;
    }

    // Treating next lines.
    for (int i = 0; i < *n; i++){
        fgets(buffer, 1025, f);
        char* strToken = strtok(buffer, "[] \n");
        int j = 0;
        while (strToken) {
            sscanf(strToken, "%lf", &(mat[i][j]));
            strToken = strtok(NULL, " ");
            j++;
        }
    }
    
    fclose(f);
    return mat;
}

int equals_matrix_2D_2D(double** A, double** B, int n){
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

    for (int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if (A[i][j] != B[i][j]){
                printf("A[%d][%d] = %f \n", i, j, A[i][j]);
                printf("B[%d][%d] = %f \n", i, j, B[i][j]);
                return 0;
            }
        }
    }

    return 1;
}

int equals_matrix_2D_1D(double** A, double* B, int n){
    /* Check if A equals B. The return values are:
    0: different
    1: equals
    */
    double** B_2D = convert_1D_to_2D(B, n);
    int res = equals_matrix_2D_2D(A, B_2D, n);
    delete_matrix_2D(&B_2D, n);
    return res;
}

int equals_matrix_file(char* filename1, char* filename2){
    /* Check if A equals B. The return values are:
    0: different
    1: equals
    */
    int n1;
    int n2;
    double** mat1 = read_matrix(filename1, &n1);
    double** mat2 = read_matrix(filename2, &n2);
    if (n1 != n2){
        delete_matrix_2D(&mat1, n1);
        delete_matrix_2D(&mat2, n2);
        return 0;
    }
    int res = equals_matrix_2D_2D(mat1, mat2, n1);
    delete_matrix_2D(&mat1, n1);
    delete_matrix_2D(&mat2, n2);

    return res;
}
