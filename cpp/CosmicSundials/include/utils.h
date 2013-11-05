#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <nvector/nvector_serial.h>
#include "psc.h"

//evaluation function that takes n_vectors as arguments
typedef int (*Fn)(N_Vector curr_guess, N_Vector residual, void *user_data);

//derivative checker. takes as input a function handle whose derivative needs
//evaluation, the point at which to evaluate the derivative, and returns the
//numerical jacobian.
cs* CheckDerivativeN(int m, int n, Fn function, N_Vector curr_guess, void* user_data);

//vector 2-norm
double Norm2(int n, double *v);

//print a vector to console
void PrintVect(int n, double *v);