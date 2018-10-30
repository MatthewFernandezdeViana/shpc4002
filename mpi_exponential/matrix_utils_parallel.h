#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h>

void convert_to_csc(
    const int *rPtr,
    const int *rcol,
    int *cPtr,
    int *crow,
    int N,
    int nnz,
    const double complex *rval,
    double complex *cval);

void generate_matrix(int *rPtr,
                     int *rCol,
                     int N,
                     int nnz,
                     double complex *rval,
                     double min,
                     double max);

void taylor_series(const double complex *val1,
                   const int *rPtr1,
                   const int *col1,
                   double complex *vector,
                   int N,
                   int nth);

double complex trace_csr(const int *rPtr1,
                         const int *col1,
                         int N,
                         const double complex *values);


void cplxd_array_gen(double complex *array, const int n, const double min, const double max);
