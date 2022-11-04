#include "lapack.hh"

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

void lapack::invert(int n,double *A) {
	int *piv = new int[n];
	double *work = new double[n];

	invert(n,A,n,piv,work,n);

	delete[] work;
	delete[] piv;
}

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
