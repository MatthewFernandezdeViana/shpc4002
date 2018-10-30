#pragma once
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>


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
                   const int nth,
                   const int loc_id,
                   const int loc_num_rows
                  );

double complex trace_csr(const int *rPtr1,
                         const int *col1,
                         int N,
                         const double complex *values);

void convert_to_csc(const int *rPtr,
                    const int *rcol,
                    int *cPtr,
                    int *crow,
                    int N,
                    int nnz,
                    const double complex *rval,
                    double complex *cval);

void distribute_mat(double complex *values,
                    int *rPtr,
                    int *colindex,
                    const int N,
                    const int loc_id,
                    const int num_procs,
                    double complex *loc_values,
                    int *loc_row_ptr,
                    int *loc_col_index,
                    int loc_num_rows
                   );

void cplxd_array_gen(double complex *array, const int n, const double min, const double max);
