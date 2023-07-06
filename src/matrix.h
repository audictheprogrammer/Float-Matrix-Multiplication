#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

// Basic operations
double** zero_matrix(int size);
double** random_matrix(int size, double p);
void delete_matrix(double*** mat, int size);
double** transpose_matrix(double** mat, int size);
void print_matrix(double** mat, int size);
void write_matrix(double** mat, int size, char* filename);
double** read_matrix(char* filename, int* size);
int equals_matrix(double** A, double** B, int size);
int equals_matrix_file(char* filename1, char* filename2);


#endif
