#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cs.h>
#include <klu.h>
#include "psc.h"
#include "utils.h"

int SolveIC(PS ps, double *u0);
int ComputeF(double *v, double *r, PS ps);
int ComputeJac(double *v, cs_di* jac, PS ps);
cs_di* InitJac(PS ps);