///////////////////////////////////////////////////////////////////////////////
// lapack.cc
//
// Class containing wrappers functions for BLAS routines
//
// Nick Derr, Nov 12, 2022
///////////////////////////////////////////////////////////////////////////////

#include "lapack.hh"


/**
 * Direct wrapper for the dgetrf_ BLAS routine, which inverts an n x n Matrix
 * in place and returns the L U factorization
 *
 * piv and work should already be allocated
 *
 * @param[in] n number of rows of square matrix
 * @param[in] A pointer to matrix elements
 * @param[in] lda lower dimension (n since assumes square)
 * @param[in] piv pointer to pivot indices
 * @param[in] work pointer to workspace array
 * @param[in] lwork length of workspace array
 */
void lapack::invert(int n,double *A,int lda,int *piv,double *work,int lwork) {

	int info;

	// LU factorize
	dgetrf_(&n,&n,A,&lda,piv,&info);

	// check for success
	if (info<0) {
		fprintf(stderr,"illegal val in %d-th arg\n",-info);
		exit(-1);
	} else if (info > 0) {
		fprintf(stderr,"U[%d,%d] = 0. A is singular\n",info,info);
		exit(-1);
	}

	// invert
	dgetri_(&n,A,&lda,piv,work,&lwork,&info);

	// check for success
	if (info<0) {
		fprintf(stderr,"illegal val in %d-th arg\n",-info);
		exit(-1);
	} else if (info > 0) {
		fprintf(stderr,"U[%d,%d] = 0. A is singular\n",info,info);
		exit(-1);
	}
}

/**
 * Inverts an n x n matrix, returning the LU factorization in place
 *
 * @param[in] n matrix size
 * @param[in] A pointer to matrix elements
 */
void lapack::invert(int n,double *A) {

	// allocate pivot and work space
	int *piv = new int[n];
	double *work = new double[n];

	// call inversion
	invert(n,A,n,piv,work,n);

	// clean up arrays
	delete[] work;
	delete[] piv;
}

/**
 * Solves the system Ax=b for x
 *
 * @param[in] n matrix size
 * @param[in] A pointer to matrix elements
 * @param[in,out] at start, b pointer to RHS.
 *   at end, pointer to solution x
 * @param[in] nrhs number of rhs columns
 */
void lapack::solve(int n,double *A,double *b,int nrhs) {
	char T = 'T';
	int lwork = 2*n*n,info;
	double *work = new double[lwork];

	dgels_(&T,&n,&n,&nrhs,A,&n,b,&n,work,&lwork,&info);

	// check for success
	if (info<0) {
		fprintf(stderr,"illegal val in %d-th arg\n",-info);
		exit(-1);
	} else if (info > 0) {
		fprintf(stderr,"U[%d,%d] = 0. A is singular\n",info,info);
		exit(-1);
	}

	delete[] work;
}

/**
 * Solves the system A x = b, given A^T
 *
 * @param[in] n matrix size
 * @param[in] AT pointer to A transpose
 * @param[in] b pointer to RHS
 * @param[out] x pointer to solution
 * @param[in] nrhs number of RHS columns
 */
void lapack::solve_T(int n,double *AT,double *b,double *x,int nrhs) {
	int *piv = new int[n];
	double *work = new double[n];
	float *swork = new float[n*(n+1)];
	int iter,info;

	dsgesv_(&n,&nrhs,AT,&n,piv,b,&n,x,&n,work,swork,&iter,&info);

	// check for success
	if (info<0) {
		fprintf(stderr,"illegal val in %d-th arg\n",-info);
		exit(-1);
	} else if (info > 0) {
		fprintf(stderr,"U[%d,%d] = 0. A is singular\n",info,info);
		exit(-1);
	}

	delete[] swork;
	delete[] work;
	delete[] piv;
}
