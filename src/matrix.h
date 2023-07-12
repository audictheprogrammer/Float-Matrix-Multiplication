#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

// Basic operations
double** zero_matrix_2D(int n);
double* zero_matrix_1D(int n);
double** random_matrix_2D(int n, double p);
double* random_matrix_1D(int n, double p);
double* convert_2D_to_1D(double** mat, int n);
double** convert_1D_to_2D(double* mat, int n);
void delete_matrix_2D(double*** mat, int n);
void delete_matrix_1D(double** mat, int n);
double** transpose_matrix(double** mat, int n);
void print_matrix_2D(double** mat, int n);
void print_matrix_1D(double* mat, int n);
void write_matrix(double** mat, int n, char* filename);
double** read_matrix(char* filename, int* n);
int equals_matrix_2D_2D(double** A, double** B, int n);
int equals_matrix_2D_1D(double** A, double* B, int n);
int equals_matrix_file(char* filename1, char* filename2);


#endif
