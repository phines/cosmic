/*
 file: KINKlu.h
 author: Ben O'Hara
 email: buohara@gmail.com
 description: Kinsol linear solver module that calls the KLU sparse direct
 linear solver. Perhaps it would be better for there to be a sparse subset of the
 KINSOL direct linear solvers of which KLU would be a specific implementation, 
 but for now I'll just attach KLU.   
 */
#ifndef KINKLU_H
#define KINLU_H

#ifdef __cplusplus
extern "C"{
#endif

//standard headers		
#include <stdio.h>
#include <stdlib.h>

//SUNDIALS headers
#include <kinsol/kinsol_impl.h>
#include <kinsol/kinsol_direct.h>
#include <sundials/sundials_direct.h>
#include <kinsol/kinsol_direct_impl.h>
#include <nvector/nvector_serial.h>

//SuiteSparse headers
#include <cs.h>
#include <klu.h>
	
/*
 * Readability replacements for accessing various pointers in the KINMem
 * memory block. These are borrowed from kin_dense.c by Radu Serban.
 */
#define lrw1           (kin_mem->kin_lrw1)
#define liw1           (kin_mem->kin_liw1)
#define func           (kin_mem->kin_func)
#define printfl        (kin_mem->kin_printfl)
#define linit          (kin_mem->kin_linit)
#define lsetup         (kin_mem->kin_lsetup)
#define lsolve         (kin_mem->kin_lsolve)
#define lfree          (kin_mem->kin_lfree)
#define lmem           (kin_mem->kin_lmem)
#define inexact_ls     (kin_mem->kin_inexact_ls)
#define uu             (kin_mem->kin_uu)
#define fval           (kin_mem->kin_fval)
#define uscale         (kin_mem->kin_uscale)
#define fscale         (kin_mem->kin_fscale)
#define sqrt_relfunc   (kin_mem->kin_sqrt_relfunc)
#define sJpnorm        (kin_mem->kin_sJpnorm)
#define sfdotJp        (kin_mem->kin_sfdotJp)
#define errfp          (kin_mem->kin_errfp)
#define infofp         (kin_mem->kin_infofp)
#define setupNonNull   (kin_mem->kin_setupNonNull)
#define vtemp1         (kin_mem->kin_vtemp1)
#define vec_tmpl       (kin_mem->kin_vtemp1)
#define vtemp2         (kin_mem->kin_vtemp2)
	
//some error messages specific to KLU
#define MSGD_KLU_JAC_GRAPH_FAIL "Sparse Jacobian structure not specified."
#define MSGD_KLU_SYMBOLIC_FAIL "KLU fill-reducing order computation failed."
#define MSGD_KLU_NUMERIC_FAIL "KLU numeric factorization failed."
#define MSGD_KLU_SOLVE_FAIL "KLU back solve failed."
	
//some numerical parameters
#define RPG_TOL 1e-5
	
	/**
	 * Attach KLU linear solver to KINSOL memory block. 
	 * 
	 * @param kin_mem Main KINSOL memory block.
	 * @param n Size (number of rows) of the problem.
	 * @param jac Initialized Jacobian.
	 * @param user_data User data to be passed to ComputeF() and ComputeJac().
	 *
	 * @return Integer status code for success/failure of attaching
	 * the solver.
	 */
	SUNDIALS_EXPORT int KINKlu(KINMem kin_memory, int num_rows);
	
	/**
	 * Provide an initialized sparse Jacobian to the solver in the form of a CXSparse sparse matrix. This function MUST be called before proceeding with a nonlinear solve. Matrices in triplet form are okay, as they'll get compressed in this function.
	 *
	 * @param jac Initialized Jacobian.
	 *
	 * @return Code indicating success of Jacobian assignment.
	 */
	SUNDIALS_EXPORT int KINKluSetJacGraph(KINMem kin_memory, cs_di *jac);
	
	/**
	 * Function prototype for the evaluation of a sparse Jacobian matrix.
	 *
	 * @param n Problem dimension.
	 * @param curr_guess Current guess for solution.
	 * @param curr_residual Residual from evaluation of current guess.
	 * @param jac A compressed sparse column matrix to hold Jacobian data.
	 * @param eval_data A pointer to any data needed to evaluate the Jacobian.
	 * @param update_fr_order Flag indicating whether the fill-reducing ordering
	 * should be recomputed. If, in the Jacobian function the user changes the
	 * matrix graph (values of row indices and column pointers), this flag should 
	 * be set to 1. 
	 * @param temp1 Temporary workspace vector for evaluating Jacobian.
	 * @param temp2 Temporary workspace vector for evaluating Jacobian.
	 *
	 * @return Status code indicating result of Jacobian evaluation.
	 */
	typedef int (*KINKluSparseJacFn)(int n, N_Vector curr_guess, 
									 N_Vector curr_residual, cs_di *jac, void *eval_data,
									 int *update_fr_order, N_Vector temp1, N_Vector temp2);
	
	/**
	 * Function prototype for matrix ordering that can be passed to KLU.
	 *
	 * @param n Number of matrix rows.
	 * @param ap Vector of column pointers (array offsets to the beginning of columns
	 * in CSC format).
	 * @param ai Vector of row indices for CSC format.
	 * @param perm A vector of permutations where perm[i] = j means that row and column
	 * j will appear as row and column i. 
	 */
	typedef int (*KINKluOrderingFn)(int n, int *ap, int *ai, 
									int *perm, klu_common *common);
	
	/**
	 A data structure for holding various KLU-related pointers like pointers
	 to Jacobian evaluation functions. This struct also holds parameters that are 
	 used by the klu_common object.
	 */
	typedef struct KINKluMemRec{
		
		///problem dimension
		unsigned int n;
		
		///pointer to a sparse jacobian evaluation function
		KINKluSparseJacFn jac_fun;
		
		/*
		 flag indicating whether the fill-reducing ordering for the jacobian
		 is up-to-date.
		 */
		unsigned int is_fill_reduced;
		
		///sparse jacobian matrix;
		cs_di *jac;
		
		///KLU common object to store various KLU linear solver parameters
		klu_common klu_comm;
		
		///KLU symbolic object used to store a fill-reducing ordering
		klu_symbolic *symbolic;
		
		///KLU numeric object to store numeric matrix factorization
		klu_numeric *numeric;
		
	} *KINKluMem;
	
	/**
	 * Data structure holding various klu-related solve parameters. These
	 * values are passed to the klu_common object stored in the KINKluMem
	 * block via the KINKluSetParams function declared below.
	 */
	typedef struct _KLUParams{
		
		///Partial pivot tolerance. Default = 0.001.
		double tol;
		
		///Which fill-reducing ordering to use. Set to 0 (AMD) by default.
		int ordering;
		
		///Whether matrix scaling should be used. Default = 2.
		int scale;
		
		///Use BTF matrix permutation. Set to 1 (true) by default.
		int btf;
		
		///Upper limit for BTF work. Default = 0 (no upper limit).
		int max_work;
		
	} *KLUParams;
	
	/**
	 * Set any non-default KLU parameters.
	 *
	 * @param kin_memory Main KINSOL memory block.
	 * @param klu_params A structure of user-defined KLU parameters.
	 */
	SUNDIALS_EXPORT int KINKluSetParams(KINMem kin_memory, KLUParams klu_params);
	
	
	/**
	 * Allows the user to specify a sparse Jacobian evaluation function for
	 * use with KLU.
	 * 
	 * @param kin_mem KINSOL memory block.
	 * @param jac_fun Pointer to a Jacobian evaluation function.
	 * 
	 * @return Status code indicating assignment result.
	 */
	SUNDIALS_EXPORT int KINKluSetSparseJacFn(void *kin_memory,
											 KINKluSparseJacFn jac_fun);
	
	/**
	 * Allow the user to specify a matrix permutation function for KLU in place
	 * of, for example, AMD. 
	 * 
	 * @param kin_mem KINSOL memory block.
	 * @param ordering_fun Pointer to a user-defined row/column permutation 
	 * function.
	 *
	 * @return Status code indicating assignment result.
	 */
	SUNDIALS_EXPORT int KINKluSetOrderingFn(void *kin_memory,
											KINKluOrderingFn ordering_fun);
	
	
#ifdef __cplusplus
}
#endif

#endif