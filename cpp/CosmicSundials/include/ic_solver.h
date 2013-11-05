/*
 file: ic_solver.h
 author: Ben O'Hara
 email: buohara@gmail.com
 description: Header file for a power system intial condition solver that uses LLNL's Sundials Kinsol nonlinear
 solver.
 */

#ifndef IC_SOLVER_H
#define IC_SOLVER_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "psc.h"
#include <cs.h>
#include <kinsol/kinsol.h>
	
/*
 depending on the presence of MPI, pull in the appropriate parallel
 or serial implementation of sundials nvector and the various
 linear solver modules.
 */
#ifdef HAVE_MPI
#include <nvector/nvector_parallel.h>
#include <kinsol/kinsol_spgmr.h>
#include <kinsol/kinsol_spbcgs.h>
#include <kinsol/sptfqmr.h>
#else
#include <nvector/nvector_serial.h>
#include <kinsol/kinsol_dense.h>
#include <kinsol/kinsol_band.h>
#include "kin_klu.h"
#endif

/**
 * Power system simulation initial condition solve.
 * 
 * @param ps Power system data.
 *
 * @return Flag indicating success/failure of initial condition solve.
 */
int SolveIC(PS ps, double *x_init, double *y_init);

					  
#ifdef __cplusplus
}
#endif
	
#endif
