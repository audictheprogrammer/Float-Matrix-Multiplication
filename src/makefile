# Targets
all: main main_test

# Headers
matrix.o: matrix.c
	gcc -c matrix.c
matrix_mul.o: matrix_mul.c
	gcc -c matrix_mul.c -fopenmp

# Mains
main: main.c matrix.o matrix_mul.o
	gcc -o main main.c matrix.o matrix_mul.o -lopenblas -lm -fopenmp
main_test: main_test.c matrix.o matrix_mul.o
	gcc -o main_test main_test.c matrix.o matrix_mul.o -lopenblas -lm -fopenmp

clean:
	rm *.o
	rm -rf main main_test
