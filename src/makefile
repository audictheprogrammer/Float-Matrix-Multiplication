# Targets
all:ain

# Headers
matrix.o: matrix.c
	gcc -c matrix.c
matrix_mul.o: matrix_mul.c
	gcc -c matrix_mul.c

# Mains
main: main.c matrix.o matrix_mul.o
	gcc -o main main.c matrix.o matrix_mul.o -lm

clean:
	rm *.o
	rm -rf main
