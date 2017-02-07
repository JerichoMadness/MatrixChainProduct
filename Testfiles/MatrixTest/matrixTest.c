#include <stdlib.h> 
#include <stdio.h>  
#include <time.h>
#include <math.h>    
#include "mkl.h"

int main() {

	double *A, *B, *C;
	int m,n,i,j,k;
	double alpha,beta;

	m = 2000;
	n = 1000;
	k = 200;

	alpha = 1.0;
	beta = 0.0;

	A = (double *)mkl_malloc( m*k*sizeof( double ), 64 );
	B = (double *)mkl_malloc( k*n*sizeof( double ), 64 );
	C = (double *)mkl_malloc( m*n*sizeof( double ), 64 );

	if (A == NULL || B == NULL || C == NULL) {
     		printf( "\n ERROR: Can't allocate memory for matrices. Aborting... \n\n");
      		mkl_free(A);
      		mkl_free(B);
      		mkl_free(C);
      		return 1;
    	}
	
	printf("Initializing Matrices...\n");

	for (i = 0; i < (m*k); i++) {
        	A[i] = (double)(i+1);
    	}	

    	for (i = 0; i < (k*n); i++) {
        	B[i] = (double)(-i-1);
    	}

    	for (i = 0; i < (m*n); i++) {
        	C[i] = 0.0;
    	}

	printf("Now computing C = A*B");
	
	for (i=0; i<100; i++) {
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                m, n, k, alpha, A, k, B, n, beta, C, n);
	}
    	printf ("\n Computations completed.\n\n");

	printf("Good job man");

	mkl_free(A);
	mkl_free(B);
	mkl_free(C);

}
