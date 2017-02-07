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
			
			if (i % 2 == 0) {		
				A[i][j] =  (((double) rand() / (double) RAND_MAX));
			} else if  (i % 2 == 1) {
				A[i][j] = -1.0* (((double) rand() / (double) RAND_MAX));
			}		
		}
	}



	
}

//Method to obtain the concrete matrix order

void MatrixChainOrder(int p[], int n, int* s){

	int i,j,k,l;
	double q;
	double m[n][n];

	double r[n+1];

	//Scale sizes of matrices to not collide with INT_MAX

	for (i=0; i<n+1; i++) {

		r[i] = (double) p[i] / (double) 1000;

	}

	//Cost is zero while multiplying one matrix
	
	for (i=0; i<n; i++)
		m[i][i] = 0;
		

	//l is chain length
	for (l=1; l<n; l++) {
		
		for (i=0; i<n-l; i++) {
		j = i+l;
		m[i][j] = INT_MAX;
			
			for (k=i; k<=j-1; k++) {
				//q = cost/scalar multiplications
				q = m[i][k] + m[k+1][j] + (p[i-1]*p[k]*p[j]);
				
				if (q < m[i][j]) {
					m[i][j] = q;
					s[n*i+j] = k;
					//printf("Yeah: s[%d][%d] %d\n", i, j, s[n*i+j]);
				}
			
			}
		}
	
	}


}




double *dgemmChainDynamic(double **A, int *s, int *p, int i, int j, int l) {

	if (i < j) {

		double *X, *Y, *Z;

		//printf("Here is the cut: %d\n", s[l*i+j]);

		X = malloc(p[i]*p[s[l*i+j]]*sizeof(double));
		
		X = dgemmChainDynamic(A,s,p,i,s[l*i+j],l);

		Y = malloc(p[s[l*i+j]+1]*p[j]*sizeof(double));

		Y = dgemmChainDynamic(A,s,p,s[l*i+j]+1,j,l);	

		Z = malloc(p[i]*p[j]*sizeof(double));

		double alpha,beta;

		alpha = 1.0; beta = 0.0;

		int m,n,k;

		m = p[i]; n = p[j]; k = p[s[l*i+j]];

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, X, k, Y, n, beta, Z, n);

		free(X);
		free(Y);
		
		//printf("Z[%d] is %1.3f\n", 4, Z[4]);

		return(Z); 

	} else if (i==j) {

		//printf("A[%d][%d] is %1.3f\n", i, 4, A[i][4]);

		return A[i];

	}

}

double* dgemmChain(double **A, int *p, int i, int j) {

	double *X, *Y, *Z;


	if (i < j) {

		X = malloc(p[i]*p[i+1]*sizeof(double));
		
		X = dgemmChain(A,p,i,i);

		Y = malloc(p[i+1]*p[j]*sizeof(double));

		Y = dgemmChain(A,p,i+1,j);	

		Z = malloc(p[i]*p[j]*sizeof(double));

		double alpha,beta;

		alpha = 1.0; beta = 0.0;

		int m,n,k;

		m = p[i]; n = p[j]; k = p[i+1];

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, X, k, Y, n, beta, Z, n);

		free(X);
		free(Y);
		
		//printf("Z[%d] is %1.3f\n", 4, Z[4]);

		return(Z); 

	} else if (i==j) {

		//printf("A[%d][%d] is %1.3f\n", i, 4, A[i][4]);

		return A[i];

	}

}


double* multiplyMatricesDynamic(double **A, int *s, int *p, int n, double *Res) {

	Res = dgemmChainDynamic(A,s,p,0,n-1,n);

	//printf("Res[%d] is %1.3f\n", 4, Res[4]);

	return Res;


}


double* multiplyMatrices(double **A, int* p, int n, double *Res) {

	Res = dgemmChain(A,p,0,n);

	return Res;

}


int main() {

	srand(time(NULL));

	double **A;
	int i,j,n;

	n = 10;

	A = (double **) malloc(n*sizeof(double*));	
	
	if (A == NULL) {
	
		printf("Error! Memory not allocated!");
		exit(0);

	}
	
	int p[n+1];

	for (i=0; i<n+1; i++) {

		p[i] = ((int)(rand() % 900)) + 100;
		//printf("p[%d]: %d\n", i,p[i]);	
	}

	initializeMatrices(p,n,A);

	if (A == NULL) {
	
		printf("Error! Memory not allocated!");
		exit(0);

	}


	printf("A[%d][%d] is %f\n", 1, 4,A[1][4]);

	int *s;

	s = malloc(n*n*sizeof(int*));

	double timerDynamic, timer;
	clock_t start, end;

	double *Res;

	Res = (double*) malloc(p[0]*p[n]*sizeof(double));

	//printf("Here is the cut: %d and the next: %d\n", s[7],s[s[7]]);

	printf("At first lets do the stuff dynamicaly\n");

	start = clock();

	MatrixChainOrder(p,n,s);

	//Res = multiplyMatricesDynamic(A,s,p,n,Res);

	end = clock();

	timerDynamic = ((double) (end - start)) / CLOCKS_PER_SEC;

	printf("The dynamic time that is needed: %1.8f\n", timerDynamic);

	printf("Now lets do the stuff normally\n");

	start = clock();

	Res = multiplyMatrices(A,p,n,Res);

	end = clock();

	timer = ((double) (end - start)) / CLOCKS_PER_SEC;

	printf("The dynamic time that is needed: %1.8f\n", timer);



	//printf(" Here: %1.3f\n",Res[4]);


	//printf("Where to split first: %d\n", s[7]);	

	/*for (i=0; i<p[0]; i++) {

		printf("[ ");
		for (j=0; j<p[n]; j++) {
			printf(" %1.10f ",Res[i*p[0]+j]);

		}
		printf(" ]\n");

	}*/

	free(A);
	free(Res);
	free(s);	
	return(0);

}
