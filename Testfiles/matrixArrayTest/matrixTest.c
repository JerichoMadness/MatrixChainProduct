#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <mkl.h>


void initializeMatrices(int p[], int n, double **A) {

	/* A is an array of pointers which point 
 	* at the start of the matrix B_x
 	* Matrices B_x are all stored in an array B
 	* Size is the length of N
 	* Start[] is the address where matrix B_x starts
 	*/

	int i,j;

	//Compute length of B

	for (i=0; i<n; i++) {

		//printf("start[%d]: %d\n", i,size);	
		A[i] = (double*) malloc(p[i]*p[i+1]*sizeof(double));
		
	}

	for (i=0; i<n; i++) {
		
		for (j=0; j<p[i]*p[i+1]; j++) {

			A[i][j] =  ((double) rand() / (double) RAND_MAX);
	
		}
	}


}

double *mult(double **A, int *p, int i, int j) {

	printf("Here is i and j: %d, %d\n", i, j);

	if (i < j) {
		
		printf("Here is i: %d\n", i);

		double *X, *Y, *Z;

		X = (double*) malloc(p[i]*p[i+1]*sizeof(double));
		
		X = mult(A,p,i,i);

		Y = (double*) malloc(p[i+1]*p[j]*sizeof(double));

		Y = mult(A,p,i+1,j);	

		Z = (double*) malloc(p[i]*p[j]*sizeof(double));

		double alpha,beta;

		alpha = 1.0; beta = 0.0;

		int m,n,k;

		m = p[i]; n = p[j]; k = p[i+1];

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, X, k, Y, n, beta, Z, n);

		//free(X);
		//free(Y);
		
		return(Z); 

	} else if (i==j) {

		printf("Here is i: %d\n",i);

		return A[i];
		
	}

}


void multiplyMatrices(double **A, int *p, int n, double *Res) {

	Res = mult(A,p,0,n-1);

}


int main() {

	srand(time(NULL));

	double **A;
	int i,n;

	n = 8;

	A = (double **) malloc(n*sizeof(double*));	
	
	if (A == NULL) {
	
		printf("Error! Memory not allocated!");
		exit(0);

	}
	
	int p[n+1];

	for (i=0; i<n+1; i++) {

		p[i] = ((int)(rand() % 100)) + 1;
		//printf("p[%d]: %d\n", i,p[i]);	
	}

	initializeMatrices(p,n,A);

	if (A == NULL) {
	
		printf("Error! Memory not allocated!");
		exit(0);

	}


	//printf("A[%d][%d] is %f\n", 1, 4,A[1][4]);

	//int *s;

	//s = malloc(n*n*sizeof(int*));

	double *Res;

	Res = malloc(p[0]*p[n]*sizeof(double));

	multiplyMatrices(A,p,n,Res);

	//printf("Where to split first: %d\n", s[7]);	

	free(A);
	free(Res);
	//free(s);	
	return(0);

}
