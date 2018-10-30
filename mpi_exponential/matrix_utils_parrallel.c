#include <stdio.h>
#include "matrix_utils_parallel.h"
#include "heap_sort.h"


/*generate an array of random numbers of different types
 * for type = 0 int, 1 double, 2 complex double.
 *
 */

int num_threads = 1;

void int_array_gen(int *array, int n, int max)
{

    int l;
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (l = 0; l < n; l++) {
        array[l] = rand() % max;
        //printf("%d\n", array[l]);
        
    }

}


/*Creates a double precision random number used for the components of
 * the complex number
 *
 */

double randfrom(int min, int max)
{
    double range = (max - min);
    double div = RAND_MAX / range;
    return min + (rand() / div);
}


/*Creates the complex array used as the values for the matrix
 *
 */

void cplxd_array_gen(double complex *array, const int n, const double min, const double max)
{

    int i;
    double real, imag;
    //double complex comp;
    
    for (i = 0; i < n; i++) {
        real = randfrom(min, max);
        imag = randfrom(min, max);
        if (rand() % 2) {
            array[i] = real + imag * I;
        } else {
            array[i] = real - imag * I;
        }
    }

}


void generate_matrix(int *rPtr,
                     int *rCol,
                     int N,
                     int nnz,
                     double complex *rval,
                     double min,
                     double max
                    )
{
    int i, j, k, num_in_col, start, stop, nplus, nminus, matrix_element = 0;
    int *col_sub_array;
    
    nplus = N + 1;
    nminus = N - 1;
    
    int_array_gen(rPtr, nplus, nnz);
    heap_sort(rPtr, nplus);
    
    rPtr[nnz] = nnz;
    
    for (i = 0; i < N; i++) {
        
        
        start = rPtr[i];
        stop = rPtr[i + 1];
        num_in_col = stop - start;
        
        //printf("%d row pointer \n", rPtr[i]);
        //printf("%d num in col\n", num_in_col);
        
        col_sub_array = calloc(num_in_col, sizeof(int));
        
        //printf("%p\n", (void*)&col_sub_array);
        
        int_array_gen(col_sub_array, num_in_col, nminus);
//
//         for (k = 0; k < num_in_col; k++){
//             //printf("%d subarray \n", col_sub_array[k]);
//             printf("%d\n", col_sub_array[k]);
//         }
        
        heap_sort(col_sub_array, num_in_col);
        
//         for (k = 0; k < num_in_col; k++){
//             printf("%d subarray \n", col_sub_array[k]);
//         }
        
        for (j = 0; j < num_in_col; j++) {
            //printf("%d column index \n", col_sub_array[j]);
            rCol[matrix_element] = col_sub_array[j];
            
            matrix_element++;
        }

        free(col_sub_array);
    }
    cplxd_array_gen(rval, nnz, min, max);
}



void convert_to_csc(
    const int *rPtr,
    const int *rcol,
    int *cPtr,
    int *crow,
    int N,
    int nnz,
    const double complex *rval,
    double complex *cval)
{
    int i, j, k, sum, temp, prev, dest;

    // count the number of values in each column
    for (i = 0; i < nnz; i++) {
        cPtr[i] = 0;
    }

    for (i = 0; i < nnz; i++) {
        cPtr[rcol[i]]++;
    }

    // Sum the counted values along column pointer
    for (j = 0, sum = 0; j < N; j++) {
        temp  = cPtr[j];
        cPtr[j] = sum;
        sum += temp;
    }
    cPtr[N] = nnz;

    // Fill out the values running through the arrays
    for (i = 0; i < N; i++) {
        for (k = rPtr[i]; k < rPtr[i + 1]; k++) {
            j  = rcol[k];
            dest = cPtr[j];

            crow[dest] = i;
            cval[dest] = rval[k];

            cPtr[j]++;
        }
    }

    for (j = 0, prev = 0; j <= N; j++) {
        temp  = cPtr[j];
        cPtr[j] = prev;
        prev    = temp;
    }

}

/*constructs the identity matrix for both csr and csc since these are the same thing
 */

void minus_mu_ident(const int *Ptr,
                    const int *index,
                    int N,
                    double complex *values,
                    double complex mu
                   )
{
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = Ptr[i]; j < Ptr[i + 1]; j++) {
            if (index[j] == i) {
                values[j] -= mu;
            }
        }
    }
}



/*Matrix norm for a csr mstrix
 * You need to re write in depth. Only a place holder since the parallel
 * one norm from:
 *
 * A block algorithm for matrix 1-norm estimation,
 * with an application to 1-norm pseudospectra
 * 2006
 *
 */


double one_norm(const int *rPtr,
                const int *rcol,
                int N,
                int nnz,
                double complex *val1)
{

    int i, k;
    int *column_pointer, *row_index;
    double complex *values_csc;
    double sum1 = 0, sum2 = 0;

    // allocate memory for the transpose used for 1 norm
    values_csc = calloc(nnz, sizeof(double complex));
    column_pointer = calloc(nnz, sizeof(int));
    row_index = calloc(N, sizeof(int));

    //fill out the transposed array
    convert_to_csc(rPtr, rcol, column_pointer, row_index, N, nnz, val1, values_csc);

    //sum across the columns
    for (i = 0; i < N; i++) {
        for (k = column_pointer[i]; k < column_pointer[i + 1]; k++) {
            sum1 += cabs(values_csc[k]);
        }
        if (sum1 > sum2) {
            sum2 = sum1;
        }
    }

    free(row_index);
    free(column_pointer);
    free(values_csc);

    return (sum2);
}


/*My algorithm for multiplication of a sparse matrix on a dense vector
 *
 */




// /*probably not essential when using Krylov methods.
//  * when reading from literature:
//  * IR is the row pointer
//  * JC is the column index
//  * Num is the value
//  */
//
// void spgemm(int *rPtrA,
//                  int *colA,
//                  int *colB,
//                  int *rPtrC,
//                  int *colC,
//                  int N,
//                 int nnzA,
//                 int nnzB,
//                  double complex *valA,
//                  double complex *valB,
//                  double complex *valC)
// {
//     int *temp_int_ptr;
//     double complex *temp_cmpx_ptr;
//     double complex *W;
//     _Bool *B;
//
//     W = calloc(N, sizeof(int));
//     B = calloc(N, sizeof(_Bool));
//
//     valC = malloc(N * sizeof(double complex));
//     rPtrC = malloc(N * sizeof(double complex));
//     valC[0] = 0;





/* Take the transpose of the second input matrix and convert it to sparse column
 */







/*Implements the trace for a sparse matrix, there's probably a faster way to do this
 * I am naive
 *
 */

double complex trace_csr(const int *rPtr1,
                         const int *col1,
                         int N,
                         const double complex *values)
{
    int i, j;
    double complex sum = 0;

    for (i = 0; i < N; i++) {
        for (j = rPtr1[i]; j < rPtr1[i + 1]; j++) {
            if (col1[j] == i) {
                sum += values[j];
                continue;
            }
            if (col1[j] > i){
                continue;
            }
        }
    }
    return (sum);

}



/* this is meant for computing the condition for
 * deciding on the scaling and squaring factor
 * however I still need to work out the
 * computation of alpha, that requires 
 * the block 1 norm estimation algorithm
 */

// int compute_p_max(m_max)
// {
//
//     float sqroot_mmax;
//     int p_low, p_high, p, pmax;
//     int i, j;
//
//     sqroot_mmax = sqrt(m_max);
//     p_low = int(floor(sqroot_mmax));
//     p_high = int(ceil(sqroot_mmax + 1));
//
//     for (p = p_low; p = p_high + 1; p++) {
//         if (p * (p - 1) <= m_max + 1) {
//             pmax = p;
//         }
//     }
//
//     return (p_max)
//
// }

unsigned long long factorial(unsigned long long n)
{
    unsigned int retval = 1;
    for (unsigned long long i = n; i > 1; --i)
        retval *= i;
    if (retval < n){
        printf("factorial size exceeded\n");
        exit(EXIT_FAILURE);        
    }
    else {
        return(retval);
    }
}

/*Taylor series evaluator to nth order
 *
 *
 */
void taylor_series(const double complex *val1,
                   const int *rPtr1,
                   const int *col1,
                   double complex *vector,                  
                   int N,
                   int nth)
{
    int i,j, fact;
    double complex *vector_next;
    
    vector_next = calloc(N, sizeof(double complex));
    
    for (j = 0; j < N; j++){
        vector_next[j] = vector[j];
    }
    
    printf("here\n");
    
    for (i = 1; i < nth + 1; i++){
        //printf("%d\n", i);
        spmv_multiply(val1, rPtr1, col1, vector, N, vector_next);
        fact = factorial(i);
        omp_set_num_threads(num_threads);
        #pragma omp parallel for
        for (j = 0; j < N; j++){
            vector[j] += vector_next[j] / fact;
        }
        
    }
    free(vector_next);

}

