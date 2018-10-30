
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>
#include "matrix_utils.h"



int main(int argc, char *argv[])
{
    int N = 10000000, NNZ = 15000000, order_approx = 10;

    MPI_Init(&argc, &argv);

    int loc_id, num_procs, loc_num_rows = false;

    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &loc_id);



    FILE *f_matrix, *f_row_pointer, *f_column_index, *f_vector, *f_result, *f_time;

//
    int min, max, i;
    int nplus, *row_pointer, *column_index;
    double complex *vector, *matrix;
    min = -1;
    max = 1;

    if (loc_id == 0) {
        column_index = calloc(NNZ, sizeof(int));
        row_pointer = calloc(N + 1, sizeof(int));
        matrix = calloc(NNZ, sizeof(double complex));
        vector = calloc(N, sizeof(double complex));

        generate_matrix(row_pointer, column_index, N, NNZ, matrix, min, max);
        cplxd_array_gen(vector, N, min, max);



        f_matrix = fopen("matrix_values.txt", "w");
        f_column_index = fopen("cindex.txt", "w");
        f_row_pointer = fopen("rpointer.txt", "w");
        f_vector = fopen("vector.txt", "w");
        f_result = fopen("result.txt", "w");
        f_time = fopen("time_omp.txt", "a");

        for (i = 0; i < NNZ; i++) {
            fprintf(f_matrix, "%f + i%f\n", creal(matrix[i]), cimag(matrix[i]));
        }
        for (i = 0; i < N + 1; i++) {
            fprintf(f_row_pointer, "%d\n", row_pointer[i]);
        }
        for (i = 0; i < NNZ; i++) {
            fprintf(f_column_index, "%d\n", column_index[i]);
        }
        for (i = 0; i < N; i++) {
            fprintf(f_vector, "%f + i%f\n", creal(vector[i]), cimag(vector[i]));
        }



        fclose(f_matrix);
        fclose(f_column_index);
        fclose(f_vector);
        fclose(f_row_pointer);
    }

    clock_t begin = clock();

    double complex *loc_values = false,  *result, *totalResult, *sendVecData;
    int *loc_row_ptr = false, *loc_col_index = false,  *loc_nnz;

    distribute_mat(matrix, row_pointer, column_index, N, loc_id, num_procs, loc_values, loc_row_ptr, loc_col_index, loc_num_rows);

    taylor_series(loc_values, loc_row_ptr, loc_col_index, vector, N, order_approx, loc_id, loc_num_rows);

    clock_t end = clock();

    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    fprintf(f_time, "%f\n", time_spent);

    for (i = 0; i < N; i++) {
        fprintf(f_result, "%f + i%f\n", creal(vector[i]), cimag(vector[i]));
    }

    fclose(f_result);


    free(vector);



}
