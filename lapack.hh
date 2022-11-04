#ifndef LAPACK_HH
#define LAPACK_HH

#include <cstdio>
#include <cstdlib>

extern "C" void dgetrf_(int *m,int *n,double *A,int *lda,int *ipv,int *info);
extern "C" void dgetri_(int *n,double *A,int *lda,int *ipiv, double *work,int *lwork,int *info);
extern "C" void dsgesv_(int *n,int *nrhs,double *A,int *lda,int *piv,double *b,int *ldb,double *x,int *ldx,double *work,float *swork,int *iter,int *info);
extern "C" void dgels_(char *T,int *m,int *n,int *nrhs,double *A,int *lda,double *b,int *ldb,double *work,int *lwork,int *info);

/**
 * A wrapper class for several LAPACK routines
 */
class lapack {

	static void invert(int n,double *A,int lda,int *piv,double *work,int lwork);

	public:

	static void invert(int n,double *A);
	static void solve_T(int n,double *A_T,double *b,double *x,int nrhs=1);
	static void solve(int n,double *A,double *b,int nrhs=1);
};

#endif
