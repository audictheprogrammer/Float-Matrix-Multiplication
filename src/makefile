# Targets
all: main main_test main_blas

# Headers
matrix.o: matrix.c
	gcc -c matrix.c
matrix_mul.o: matrix_mul.c
	gcc -c matrix_mul.c

# Mains
main: main.c matrix.o matrix_mul.o
	gcc -o main main.c matrix.o matrix_mul.o -lm
main_test: main_test.c matrix.o matrix_mul.o
	gcc -o main_test main_test.c matrix.o matrix_mul.o -lm
main_blas: main_blas.c
	gcc -o main_blas main_blas.c -lopenblas

clean:
	rm *.o
	rm -rf main main_test main_blas
