/*
 file: KINKlu.c
 author: Ben O'Hara
 email: buohara@gmail.com
 description: Implementation of a KINSOL interface that uses KLU
 for sparse linear solves.
 */

#include "kin_klu.h"

/**
 * Initialize and allocate memory for KLU.
 *
 * @param kin_mem KINSOL memory block.
 * 
 * @return Integer status code for result of initialization.
 */
static int InitKINKlu(KINMem kin_memory);
	
/**
 * Setup KLU for a linear solve, namely, do symbolic analysis
 * numeric (LU) factorization.
 *
 * @param kin_mem KINSOL memory block.
 *
 * @return Integer status code for result of setup.
 */
static int SetupKINKlu(KINMem kin_memory);
	
/**
 * Perform the sparse linear solve using KLU.
 *
 * @param kin_mem KINSOL memory block.
 * @param x Solution vector to the linear system Ax=b.
 * @param b Right-hand side input of Ax=b.
 * @param res_norm The norm of the residual from the previous guess for x.
 *
 * @return Integer status code for result of solve.
 */
static int SolveKlu(KINMem kin_memory, N_Vector x, N_Vector b, realtype *res_norm);
	
/**
 * Free up memory and resources used by KLU.
 *
 * @param kin_mem KINSOL memory block.
 */
static void FreeKINKlu(KINMem kin_memory);

/*
 exported function that attaches the KLU linear solver in a
 main program. this initializes the kinsol memory block as well as a memory
 block for klu-related parameters and functions.
 */
int KINKlu(KINMem kin_memory, int num_rows){
	
	//kinsol and kinsol sparse memory blocks
	KINKluMem kin_klu_mem;
	
	/* 
	 Return immediately if kinmem is NULL 
	 */
	if (!kin_memory) {
		KINProcessError(NULL, KINDLS_MEM_NULL, "KINKLU", "KINKlu", MSGD_KINMEM_NULL);
		return(KINDLS_MEM_NULL);
	}
	
	/*
	 Free the Kinsol linear solver memory block if it hasn't already been freed.
	 */
	if(kin_memory->kin_lfree){ 
		kin_memory->kin_lfree(kin_memory);
	}
	
	/*
	 attach the four KLU solver interface functions to the appropriate
	 function pointers in the KINSOL memory block
	 */
	kin_memory->kin_linit = InitKINKlu;
	kin_memory->kin_lsetup = SetupKINKlu;
	kin_memory->kin_lsolve = SolveKlu;
	kin_memory->kin_lfree = FreeKINKlu;
	
	/*
	 Next, we allocate room for a KINKlu memory block and attach the memory block 
	 to the linear solver pointer of the main kinsol memory block.
	 */
	kin_klu_mem = malloc(sizeof(struct KINKluMemRec));
	
	//check for memory error
	if(!kin_klu_mem) return 1;
	
	//set the problem size in the kinklu memory block
	kin_klu_mem->n = num_rows;
	
	//attach the klu memory block to the main solver block
	kin_memory->kin_lmem = kin_klu_mem;
	
	//return
	return 0;
};

/*
 Initialize a sparse Jacobian for the solver.
 */ 
int KINKluSetJacGraph(KINMem kin_memory, cs_di *jac){

	//check inputs
	if(!kin_memory||!jac) return 1;

	//grab the linear solver memory block
	KINKluMem kin_klu_mem=(KINKluMem)kin_memory->kin_lmem;
	if(!kin_klu_mem) return 1;

	//grab needed kinklu objects
	cs_di *kin_klu_jac=kin_klu_mem->jac;
	
	//check whether the jacobian data has already been allocated. if so, free it.
	if(kin_klu_jac){
		cs_di_spfree(kin_klu_jac);
		kin_klu_jac=NULL;
	}

	//check if the input jacobian is in triplet form. compress if so.
	if(CS_TRIPLET(jac)){
		kin_klu_jac=cs_di_compress(jac);
		if(!kin_klu_jac) return 1;
		cs_di_spfree(jac);
		jac = kin_klu_jac;
		return 0;
	}
	
	//assign the input jacobian to the kinklu jacobian
	kin_klu_mem->jac=jac;
	
	//return success
	return 0;
};

/*
 This will prep KLU for subsequent linear solves by computing an initial fill-reducing ordering and numeric factorization (if values are available) of the user-supplied Jacobian.
 */
int InitKINKlu(KINMem kin_memory){
	
	//check inputs
	if(!kin_memory) return 1;
	
	//grab the kinklu block
	KINKluMem kin_klu_mem=(KINKluMem)kin_memory->kin_lmem;
	if(!kin_klu_mem) return 1;
	
	kin_memory->kin_inexact_ls=FALSE;
	kin_memory->kin_setupNonNull=TRUE;
	
	//grab klu objects from kinklu memory
	int n=kin_klu_mem->n, ok;
	cs_di *jac=kin_klu_mem->jac;
	klu_symbolic *symb=kin_klu_mem->symbolic;
	klu_numeric *numeric=kin_klu_mem->numeric;
	klu_common *comm=&(kin_klu_mem->klu_comm);
	
	//if theres no jacobian to try factoring, we're done
	if(!jac) return 1;
	
	//check for, and delete if necessary, existing klu symbolic and numeric objects.
	if(symb){klu_free_symbolic(&symb, comm); symb=NULL; kin_klu_mem->symbolic=NULL;}
	if(numeric){klu_free_numeric(&numeric, comm); numeric=NULL; kin_klu_mem->numeric=NULL;}
	ok=klu_defaults(comm);
	
	//attempt the symbolic factorization.
	symb=klu_analyze(n, jac->p, jac->i, comm);
	
	//if there's an error doing the factorization, abort. don't delete the jacobian passed in, but null the kinklu jac pointer
	if(!symb) return 1;
	
	//otherwise, assign the kinklu symbolic pointer and indicate that the jacobian has been fill-reduced
	kin_klu_mem->symbolic=symb;
	kin_klu_mem->is_fill_reduced=1;
	
	//do numeric factor if values are available
	/*
	if(jac->x){
		//try factorization. on failure, free the previous symbolic object and return.
		numeric = klu_factor(jac->p, jac->i, jac->x, symb, comm);
		if(!numeric){
			klu_free_symbolic(&symb, comm);
			kin_klu_mem->symbolic = NULL;
			return 1;
		}
		//otherwise, assign the kinklu numeric pointer and return.
		kin_klu_mem->numeric = numeric;
	}*/
	return 0;
};

/*
 Setup KLU for a linear solve. This function factors the Jacobian matrix
 before handing off the factors for a back solve. Optionally, this function also
 computes a new fill-reducing ordering (using KLU) in the case that the matrix 
 graph has been updated.
 */
int SetupKINKlu(KINMem kin_memory){
	
	//get the KINKlu memory block
	KINKluMem kin_klu_mem=(KINKluMem)kin_memory->kin_lmem;
	if(!kin_klu_mem) return 1;
	
	//grab appropriate klu objects
	cs_di *jac=kin_klu_mem->jac;
	klu_symbolic *symb=kin_klu_mem->symbolic;
	klu_numeric *numeric=kin_klu_mem->numeric;
	klu_common *comm=&(kin_klu_mem->klu_comm);
	int n=kin_klu_mem->n, update_fr_order=0;
	
	//call the jacobian evaluation function
	kin_klu_mem->jac_fun(n, kin_memory->kin_uu, kin_memory->kin_fval, jac, kin_memory->kin_user_data, &update_fr_order, kin_memory->kin_vtemp1, kin_memory->kin_vtemp2);
	
	/*
	 if a new fill-reducing ordering has been requested, or if the graph and values have been specified but no ordering has been computed yet, perform the computation
	 */
	if(update_fr_order){
	
		//if a symbolic object already exists, free it
		if(symb){
			klu_free_symbolic(&symb, comm);
			kin_klu_mem->symbolic=NULL;
		}
		
		//perform the fill-reducing ordering
		symb=klu_analyze(n, jac->p, jac->i, comm);
		if(!symb) return 1;
		kin_klu_mem->symbolic=symb;
		
		/*
		now we need to perform a numeric factorization. first, free an existing
		numeric factorization if there is one.
		*/
		if(numeric){
			klu_free_numeric(&numeric, comm);
			kin_klu_mem->numeric=NULL;
		}
		
		//perform the factorization
		numeric=klu_factor(jac->p, jac->i, jac->x, symb, comm);
		
		/*
		 check if the factorization was successful and return if not
		 */
		if(!numeric){
			klu_free_symbolic(&symb, comm);
			kin_klu_mem->symbolic=NULL;
			return 1;
		}
		kin_klu_mem->numeric=numeric;
		
		//otherwise, the factorization was a success and we can return
		return(KINDLS_SUCCESS);		
	}
	
	/*
	 if a new fill-reducing ordering is not necessary, we can proceed with factorization. first, check if a numeric factorization exists. if not, compute it
	 */
	if(!numeric)
	{
		//perform the factorization
		numeric=klu_factor(jac->p, jac->i, jac->x, symb, comm);
		
		/*
		 check if the factorization was successful and return if not
		 */
		if(!numeric) return 1;
		kin_klu_mem->numeric=numeric;
		return(KINDLS_SUCCESS);
	}
	
	/*
	 if a symbolic and numeric factorization already exist, try a refactor using the old numeric factorization. this is much faster than a full numeric factorization and requires no new memory
	*/
	klu_refactor(jac->p, jac->i, jac->x, symb, numeric, comm);
		
		
#ifdef _VERBOSE
		
	/*
	 check the pivot growth factor. i confess that i dont understand what this
	 factor means, but a small value is supposed to indicate numerical
	 instability for testing i'm going to compute it and do a full numerical
	 factorization if it's too small.
	 */
	klu_rgrowth(jac->p, jac->i, jac->x, symb, numeric, comm);
		
	/*
	 print the growth factor to the console for testing
	 */
	printf("Reciprocal pivot growth after refactor: %1.5e\n\n", kin_klu_mem->klu_comm.rgrowth);
		
#endif
	
	//return
	return(KINDLS_SUCCESS);
};

/*
 Perform linear solve with KLU. This process is pretty quick once the
 fill-reducing orderings and numeric factorizations are performed.
 */
int SolveKlu(KINMem kin_memory, N_Vector x, N_Vector b, realtype *res_norm){
		
	//get the klu memory block
	KINKluMem kin_klu_mem = (KINKluMem)kin_memory->kin_lmem;
	if(!kin_klu_mem) return 1;
	
	//grab the right-hand side data from the nvector passed in
	double *b_data = NV_DATA_S(b);
	
	int i,n=kin_klu_mem->n;
	klu_symbolic *symb=kin_klu_mem->symbolic;
	klu_numeric *numeric=kin_klu_mem->numeric;
	klu_common *comm=&(kin_klu_mem->klu_comm);
	
	SetupKINKlu(kin_memory);
	
	//perform the linear solve
	klu_solve(symb, numeric, n, 1, b_data, comm);
	NV_DATA_S(x)=b_data;

	//return
	return 0;
};

/*
 Free memory used by KLU
 */
void FreeKINKlu(KINMem kin_memory){
	
	/*
	do memory deallocation from the bottom up. first, deallocate memory
	used by the KLU memory block.
	*/
	
	//grab the KLU memory block
	KINKluMem kin_klu_mem = (KINKluMem)kin_memory->kin_lmem;
	if(!kin_klu_mem) return;
	
	//free the jacobian matrix
	if(kin_klu_mem->jac) cs_di_spfree(kin_klu_mem->jac);
	
	//free the klu objects
	if(kin_klu_mem->symbolic) klu_free_symbolic(&(kin_klu_mem->symbolic), &(kin_klu_mem->klu_comm));
	
	if(kin_klu_mem->numeric) klu_free_numeric(&(kin_klu_mem->numeric), &(kin_klu_mem->klu_comm));
	
	//now release the klu memory block
	free(kin_klu_mem);
	kin_memory->kin_lmem = NULL;
	
	//return
	return;
};

/*
 * Allows the user to specify a sparse Jacobian evaluation function for
 * use with KLU.
 */
int KINKluSetSparseJacFn(void *kin_memory, KINKluSparseJacFn jac_fun){
	
	//create KINSOL and KINKlu memory blocks.
	KINMem kin_mem;
	KINKluMem kin_klu_mem;
	
	//check that the KINSOL memory block is not null.
	if (!kin_memory) {
		KINProcessError(NULL, KINDLS_MEM_NULL, "KINKLU", "KINKluSetSparseJacFn", MSGD_KINMEM_NULL);
		return(KINDLS_MEM_NULL);
	}
	
	//grab it if okay
	kin_mem=(KINMem)kin_memory;
	
	//check if the klu memory block is null
	if(!kin_mem->kin_lmem){
		KINProcessError(NULL, KINDLS_MEM_NULL, "KINKLU", "KINKluSetSparseJacFn", MSGD_KINMEM_NULL);
		return(KINDLS_MEM_NULL);
	}
	
	//if not, grab it
	kin_klu_mem=(KINKluMem)kin_mem->kin_lmem;
	
	//now attach the jacobian function to the klu memory block
	kin_klu_mem->jac_fun=jac_fun;
	
	//return
	return(KINDLS_SUCCESS);
}

/*
 * Allow the user to specify a matrix permutation function for KLU in place
 * of, for example, AMD. 
 */
int KINKluSetOrderingFn(void *kin_memory, KINKluOrderingFn ordering_fun){
	
	//create KINSOL and KINKlu memory blocks.
	KINMem kin_mem;
	KINKluMem kin_klu_mem;
	
	//check that the KINSOL memory block is not null.
	if(!kin_memory){
		KINProcessError(NULL, KINDLS_MEM_NULL, "KINKLU", "KINKluSetOrderingFn", MSGD_KINMEM_NULL);
		return(KINDLS_MEM_NULL);
	}
	
	//cast it if it's not
	kin_mem = (KINMem) kin_memory;
	
	//check if the klu memory block is null
	if(!kin_mem->kin_lmem){
		KINProcessError(NULL, KINDLS_MEM_NULL, "KINKLU", "KINKluSetOrderingJacFn", MSGD_KINMEM_NULL);
		return(KINDLS_MEM_NULL);
	}
	
	//if not, grab it
	kin_klu_mem = (KINKluMem)kin_mem->kin_lmem;
	
	//set the klu memory block's ordering function to the one passed in
	kin_klu_mem->klu_comm.user_order = ordering_fun;
	
	//return with triumphant success
	return(KINDLS_SUCCESS);
}

/*
 * Since we are using sparse matrices, we need to know the structure of the
 * nonzeros before we can peform any solves. If, in advance, the user wishes
 * to specify the nonzero structure of the matrix to be used, he/she can do
 * so here (if, for example, the graph remains constant and only the values
 * change). The values in the arrays passed in will be copied, so the user is
 * free to delete the arrays from user code after calling this function.
 */
 
 /*
int KINKluSetGraph(void *kin_memory, int nnz, 
				   int *col_pointers, int *row_indices){
	
	//create KINSOL and KINKlu memory blocks.
	KINMem kin_mem;
	KINKluMem kin_klu_mem;
	
	//check that the KINSOL memory block is not null.
	if (kin_memory == NULL) {
		KINProcessError(NULL, KINDLS_MEM_NULL, "KINKLU", 
						"KINKluSetGraph", MSGD_KINMEM_NULL);
		return(KINDLS_MEM_NULL);
	}
	
	//cast it if it's not
	kin_mem = (KINMem) kin_memory;
	
	//check if the klu memory block is null
	if(kin_mem->kin_lmem == NULL){
		KINProcessError(NULL, KINDLS_MEM_NULL, "KINKLU", 
						"KINKluSetGrapg", MSGD_KINMEM_NULL);
		return(KINDLS_MEM_NULL);
	}
	
	//if not, grab it
	kin_klu_mem = (KINKluMem)kin_mem->kin_lmem;
	
	//get the problem size
	int num_cols = kin_klu_mem->n;
	
	//
	 //the number of column pointers will be num_cols+1, so allocate
	//this much memory and copy the column pointers.
	//
	int i;
	if(kin_klu_mem->jac->col_pointers != NULL)
		free(kin_klu_mem->jac->col_pointers);
	kin_klu_mem->jac->col_pointers = calloc(num_cols + 1, sizeof(int));
	for(i = 0; i < num_cols + 1; i++)
		kin_klu_mem->jac->col_pointers[i] = col_pointers[i];
	
	//
	 //Similarly, the number of row indices will be nnz, so allocate this
	 //memory and copy.
	 //
	if(kin_klu_mem->jac->row_indices != NULL)
		free(kin_klu_mem->jac->row_indices);
	kin_klu_mem->jac->row_indices = calloc(nnz, sizeof(int));
	for(i = 0; i < nnz; i++)
		kin_klu_mem->jac->row_indices[i] = row_indices[i];
	
	//
	 //Now that the matrix graph is specified, we can compute a fill-reducing
	 //ordering, so invoke klu_analyze.
	 //
	kin_klu_mem->symbolic = klu_analyze(num_cols, kin_klu_mem->jac->col_pointers,
										kin_klu_mem->jac->row_indices, &(kin_klu_mem->klu_comm));
	
	//
	 //check if the computation was successful
	 //
	if(kin_klu_mem->symbolic == NULL)
	{
		//process the error and return
		KINProcessError(kin_memory, KINDLS_JACFUNC_UNRECVR, "KINKLU", "KINKluSetGraph", 
						"KLU symbolic factorization failed.");
		return(KINDLS_JACFUNC_UNRECVR);
		
	}
	
	//
	 //Now that we have an initial fill-reducing ordering, we can set the
	 //is_fill_reduced flag in klu memory to true.
	 //
	kin_klu_mem.is_fill_reduced = 1;
	
	//return
	return(KINDLS_SUCCESS);
}
*/

/*
 Set any non-default KLU parameters.
*/
int KINKluSetParams(KINMem kin_memory, KLUParams klu_params){
	
	//grab the klu memory block
	KINKluMem kin_klu_mem = (KINKluMem)kin_memory->kin_lmem;
	
	//and then fill the klu_common object with the values specified.
	kin_klu_mem->klu_comm.tol = klu_params->tol;
	kin_klu_mem->klu_comm.ordering = klu_params->ordering;
	kin_klu_mem->klu_comm.scale = klu_params->scale;
	kin_klu_mem->klu_comm.btf = klu_params->btf;
	kin_klu_mem->klu_comm.maxwork = klu_params->max_work;
	
	//and return success
	return(KIN_SUCCESS);
}