#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

// Basic operations
float** zero_matrix(int size);
float** random_matrix(int size, int p);
void delete_matrix(float*** mat, int size);
float** transpose_matrix(float** mat, int size);
void print_matrix(float** mat, int size);
void write_matrix(u_int64_t** mat, int size, char* filename);
u_int64_t** read_matrix(char* filename, int* size);
int equals_matrix(u_int64_t** A, u_int64_t** B, int size);
int equals_matrix_file(char* filename1, char* filename2);


#endif
