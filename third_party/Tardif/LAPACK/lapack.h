
/* extern void dgemm_(char*,char*, int*, int*,int*, double*, double*, int*, double*, int*, double*, double*, int*); */

extern  int dgemm_(char *transa, char *transb, int *m, int *n,
		   int *k, double *alpha, double *a, int *lda,
		   double *b, int *ldb, double *beta, double *c, int *ldc);

#include "../CMINPACK/f2c.h"

/* extern /\* Subroutine *\/ int dgemm_(char *, char *, integer *, integer *, */
/* 				   integer *, doublereal *, doublereal *, integer *, doublereal *, */
/* 				   integer *, doublereal *, doublereal *, integer *); */

extern int dgemv_(char *trans, int *m, int *n, double *alpha,
		  double *a, int *lda, double *x, int *incx,
		  double *beta, double *y, int *incy);


/* extern int dsyev_(char *jobz, char *uplo, int *n, double *a, */
/* 		  int *lda, double *w, double *work, int *lwork, */
/* 		  int *info); */


/* extern int dsyev_(char *N, char* L,int* order, double* mJmJt, int* LDA, double* eigVal, double* work, int* lwork, int* info); */
