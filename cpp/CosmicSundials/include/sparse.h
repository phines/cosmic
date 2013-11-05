/*
 file: sparse.h
 author: Ben O'Hara
 email: buohara@gmail.com
 description: Implementation of various sparse data structures and
 operations.
 */

#ifndef SPARSE_H
#define SPARSE_H

#ifdef __cplusplus
extern "C"{
#endif
	
#include <stdio.h>
#include <stdlib.h>

	///Admittance matrix structure.
	typedef struct{
	
		///Flag indicating whether matrix entries are complex.
		int is_complex; 
	
		///Matrix dimensions.
		unsigned int num_rows, num_cols;
	
		///Column indices and row pointers.
		int *row_indices;
		int *col_pointers;
	
		///Values.
		double *vals;	
	}CSCMatrix, *CSCMat;

	///A general matrix of length n.
	typedef struct{
	
		///Flag indicating whether vector entries are complex.
		int is_complex;
	
		///Vector length
		unsigned int n;
	
		///Values
		double *vals;
	}Vector, *Vec;
	
	/**
	* Sparse matrix-vector multiplication.
	*
	* @param matrix A general matrix in compressed sparse row format.
	* @param vec General vector of length n.
	*/
	void CSCMV(CSCMat matrix, Vec in_vec, Vec out_vec);

	/**
	 * Free resources used by a CSC matrix.
	 *
	 * @param matrix Pointer to a CSC matrix to be freed.
	 */
	void FreeCSCMatrix(CSCMat matrix);
	
#ifdef __cplusplus
}
#endif

#endif