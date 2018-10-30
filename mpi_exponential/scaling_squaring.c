#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <time.h>

#include "matrix_utils.h"

bool condition(double complex A_norm, int m_max, int columns_norm_est=5)
{

    int a;
    double b;
    int pmax = compute_p_max(m_max);


    a = 2 * columns_norm_est * p_max * (p_max + 3);
    b = theta[m_max] / double(m_max);

    return (A_norm <= a * b)
}

void m_s_frag3_1(int *rPtr1,
                 int *col1,
                 int N,
                 int nnz,
                 double complex *values,
                 double tol,
                 int m_fin,
                 int s_fin,
                 double complex norm,
                 int m_max = 55
                )
{
    int m, p, m_fin, s_fin;

    // these thetas had to be computed since the papers don't give them
    double theta[55] = {0.000000000000000000001,
                        0.0000000001,
                        0.000013863478661191213,
                        0.0003397168839976962,
                        0.002400876357887274,
                        0.009065656407595102,
                        0.023844555325002736,
                        0.049912288711153226,
                        0.08957760203223342,
                        0.14418297616143777,
                        0.21423580684517107,
                        0.2996158913811581,
                        0.3997775336316795,
                        0.5139146936124294,
                        0.6410835233041199,
                        0.7802874256626575,
                        0.9305328460786568,
                        1.0908637192900361,
                        1.2603810606426387,
                        1.438252596804337,
                        1.6237159502358214,
                        1.8160778162150857,
                        2.014710780944616,
                        2.21904886936509,
                        2.4285825244428265,
                        2.6428534574594353,
                        2.861449633934264,
                        3.084000544989162,
                        3.310172839890271,
                        3.5396663487436895,
                        3.772210495681751,
                        4.00756108611804,
                        4.245497442579696,
                        4.485819859447368,
                        4.728347345793539,
                        4.972915626191982,
                        5.219375371084059,
                        5.467590630524544,
                        5.717437447572013,
                        5.968802630041848,
                        6.221582661689891,
                        6.475682736079984,
                        6.731015898381024,
                        6.98750228213063,
                        7.245068429597952,
                        7.503646685788864,
                        7.763174657377988,
                        8.02359472893998,
                        8.284853629803916,
                        8.546902045684934,
                        8.809694269971322,
                        9.073187890176143,
                        9.337343505612015,
                        9.602124472826556,
                        9.8674966757534
                       }

//     if (condition(norm, m_max) == 1) {
//         m_fin = 0; //TODO set the initial values for m and s
//         for (m = 0; m < 55; m++) {
//             s = int(ceil(norm / theta[m]));
//             if (m * s < m_fin * s_fin) {
//                 fin_m = m
//                         fin_s = s
//             }
//         }
//     } else {
//         for (p = 2; p = compute_p_max(m_max) + 1; p++) {
//             for (m = p * (p - 1) - 1) {
//
//             }
//         }
//     }



}


/*Here lies the matrix exponential algorithm, balancing is not done :( this is
 * so sad, Alexa play Despacito. It seems like the papers
 * you cheated and used the block 1 norm in python, this turns out to be hardish
 * 
 */

void exp_matrix_vector(double complex *val1,
                       int *rPtr1,
                       int *col1,
                       double complex *vector,
                       int N,
                       int nnz,
                       double complex *result)
{
    //this designates the algorithm tolerance
    // long int tol = pow((long int) 2, -24); //uncomment for single precision
    long int tol = pow((long int) 2, -53);    //comment out for single
    double complex mu, A1norm, m, s;

    mu = trace_csr(rptr1, col1, N, nnz) / N;
    minus_mu_ident(rPtr1, col1, N, val1, mu);
    A1norm = one_norm(rPtr1, col1, N, nnz, val1);

    if (0 == A1norm) {
        m = 0;
        s = 1;
    } else {
        parameters(m, s, tol);
    }
}


/*Main algorithm to implement Arnoldi iterations to get eigenvalues estimates.
 *this imports the matrices and runs over them creating the krylov subspace
 * however you still need to exponentiate the upper Hessenberg matrix at the end
 * (╯°□°）╯︵ ┻━┻
 * TODO pull code from your python implementation
 * 
 * 
 *
 *
 */

int main()
{
    srand((unsigned) time(NULL));

    int An, Annz, Aindex, Aptr;
    double complex Avalues, Bvalues, Bvec, Fvec;


    //TODO read off values of An and Annz from file







}
