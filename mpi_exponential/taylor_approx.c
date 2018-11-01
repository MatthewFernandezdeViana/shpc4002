
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
    

    MPI_Init(&argc, &argv);
    int N = atoi(argv[1]), NNZ = atoi(argv[2]), order_approx = atoi(argv[3]); //define the matrix size

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

    column_index = calloc(NNZ, sizeof(int));
    row_pointer = calloc(N + 1, sizeof(int));
    matrix = calloc(NNZ, sizeof(double complex));
    vector = calloc(N, sizeof(double complex));

    printf("not cooking yet\n");

    if (loc_id == 0) {

        double tzero, t1, t2, t3, t4, t5, tfinal;


        generate_matrix(row_pointer, column_index, N, NNZ, matrix, min, max); //make the matrices to be acted upon
        cplxd_array_gen(vector, N, min, max);

        printf("here we are friend\n");
        
        // the following code id for saving and writing the matrix and original vector hwoever 
        // it is commented out so that it doesn't write when testing or you don't desire the full thing



//         f_matrix = fopen("matrix_values.txt", "w");
//         f_column_index = fopen("cindex.txt", "w");
//         f_row_pointer = fopen("rpointer.txt", "w");
//         f_vector = fopen("vector.txt", "w");
//         f_result = fopen("result.txt", "w");
//         f_time = fopen("time_omp.txt", "a");
//
//         printf("Files opened\n");
//
//         for (i = 0; i < NNZ; i++) {
//             fprintf(f_matrix, "%f + i%f\n", creal(matrix[i]), cimag(matrix[i]));
//         }
//         for (i = 0; i < N + 1; i++) {
//             fprintf(f_row_pointer, "%d\n", row_pointer[i]);
//         }
//         for (i = 0; i < NNZ; i++) {
//             fprintf(f_column_index, "%d\n", column_index[i]);
//         }
//         for (i = 0; i < N; i++) {
//             fprintf(f_vector, "%f + i%f\n", creal(vector[i]), cimag(vector[i]));
//         }
//
//
//
//         fclose(f_matrix);
//         fclose(f_column_index);
//         fclose(f_vector);
//         fclose(f_row_pointer);
//         printf("matrix written\n");
    }

    if (loc_id == 0) {
        tzero = MPI_Wtime();
    }

    double complex *loc_values = false,  *result, *totalResult, *sendVecData;
    int *loc_row_ptr = false, *loc_col_index = false,  *loc_nnz;
    
    
    /*Send out all the matrix values and distribute them by sending a number of rows depending on the process number     * 
     * 
     */
    if (num_procs != 1){
    distribute_mat(matrix, row_pointer, column_index, N, loc_id, num_procs, loc_values, loc_row_ptr, loc_col_index, loc_num_rows);
    }

    if (loc_id == 0) {
        t1 = MPI_Wtime();
    }

    //Perform the Taylor series expansion to the order of size specified 
    taylor_series(loc_values, loc_row_ptr, loc_col_index, vector, N, order_approx, loc_id, loc_num_rows);

    if (loc_id == 0) {
        tfinal = MPI_Wtime();
    }

    //clock_t end = clock();

    if (loc_id == 0) {

        double time_spent = (double)(tzero - tfinal) / CLOCKS_PER_SEC;
        fprintf(f_time, "%f\n", time_spent);

        for (i = 0; i < N; i++) {
            fprintf(f_result, "%f + i%f\n", creal(vector[i]), cimag(vector[i]));
        }

        fclose(f_result);
    }


    free(vector);



}
