/*
 file: sparse.c
 author: Ben O'Hara
 email: buohara@gmail.com
 description: Implementation of some basic sparse linear algebra routines.
 */

#include "sparse.h"

//sparse matrix-vector multiplication
void CSCMV(CSCMat matrix, Vec in_vec, Vec out_vec){
	
	//return
	return;
}

//free up data used by this matrix.
void FreeCSCMatrix(CSCMat matrix){
	
	//just free the pointers in this matrix
	free(matrix->row_indices);
	free(matrix->col_pointers);
	free(matrix->vals);
	
	//return
	return;
}