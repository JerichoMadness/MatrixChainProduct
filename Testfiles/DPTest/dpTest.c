#include <stdlib.h> 
#include <stdio.h>  
#include <time.h>
#include <math.h>
#include <limits.h>    
#include <mkl.h>


double Mult(A,s,i,j) {

	if (i<j) {

		double X = Mult(A,s,i,s[i,j]);
		double Y = Mult(A,s,s[i+j]+1,j);

		return X*Y;
	
	}

	else return A[i];

}

//Method to obtain the concrete matrix order

int MatrixChainOrder(int p[], int n) {

	int m[n][n];
	int s[n][n]
	int i,j,k,l,q;
	

	//Cost is zero while multiplying one matrix
	
	for (i=1; i<n; i++)
		m[i][i] = 0;
		

	//l is chain length
	for (l=2; l<n; l++) {
		
		for (i=1; i<n-l+1; i++) {
		j = i+l-1;
		m[i][j] = INT_MAX;
			
			for (k=i; k<=j-1; k++) {
				//q = cost/scalar multiplications
				q = m[i][k] + m[k+1][j] + p[i-1]*p[k]*p[j];
				if (q < m[i][j])
					m[i][j] = q;
			
			}
		}
	
	}

	return m[1][n-1];

}


int main() {

	int arr[] = {40, 20, 30, 10, 30};
	int size = sizeof(arr)/sizeof(arr[0]);
	
	printf("Minimum num of multiplications is %d", MatrixChainOrder(arr,size));

	getchar();
	
	return 0;

}
