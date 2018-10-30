# The compiler: gcc in this case, will change when moving to MPI
CC = mpicc

# Compiler flags:
CFLAGS = -std=c11 -O0 -fopenmp -Wall -g
LINKERS = -lm

TARGET = shpc2
SRCS_PAR = taylor_approx.c matrix_utils_parrallel.c heap_sort.c
SRCS_SER = taylor_approx.c matrix_utils.c heap_sort.c
HEADERS = matrix_utils.h heap_sort.h
HEADERS_par = matrix_utils_parallel.h heap_sort.h

build: $(SRCS)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRCS_SER) $(HEADERS) $(LINKERS)
	#$(CC) $(CFLAGS) -o $(TARGET)_parallel $(SRCS_PAR) $(HEADERS_par) $(LINKERS)

clean:
	rm -f *.o
	rm $(TARGET)
	rm $(TARGET)_parallel

rebuild: clean build
