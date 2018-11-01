 #!/bin/bash
 module swap PrgEnv-cray PrgEnv-gnu
 make
read -p 'Please put in matrix size:' N
read -p 'Please put in number of non zero entries:' nnz
read -p 'Please enter order of approximation, preferably less than 13' order_approximation

