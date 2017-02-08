# MatrixChainProduct
Final project for HPMC

Task:

** 1) Matrix chain product
- X := M_1 M_2 ... M_n, what is the optimal parenthesization?
- classical solution: minimum flops
  best known algorithm: O(n log n)   <- forget about it
  standard dynamic programming algorithm: O(n^3)  
- modify the standard O(n^3) algorithm to find the solution for the best execution time
- the matrix products are performed through BLAS calls
- solve both for single core and multicore
- run comparisons
- find examples where the solution for best time and for minimum flops lead to a large difference in execution time

Structure:

Testfiles:
	Random files for testing purposes. Not relevant

Code: 
	All the main code files

Code/BatchFiles:
	Scripts for the cluster

Code/BatchFiles/RandomTests:
	Output files for random matrix size entries

Code/BatchFiles/ExplicitTests:
	Output files for specific matrix entries
