/* SNAPPER core functions
   Copyright (C) 1991-2020 Free Software Foundation, Inc.
   This file is part of SNAPPER.

   SNAPPER is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   SNAPPER is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <https://www.gnu.org/licenses/>. 
     
    Authors: Marnus Stoltz (marnus@gbsassociates.co.za)
             David Bryant (david.bryant@otago.ac.nz)
 */

#include "snapper_core.h"

// Number of poles used in Caratheodory-Fejer procedure to approximate matrix exponent
#if (CF_DEG == 12)
// Constant used in Caratheodory-Fejer procedure to approximate matrix exponent
const double complex c[] = {0.000818433612497 + 0.000581353207069*I,
	                       -0.068571505514864 - 0.038419074245887*I,
	                        1.319411815998137 + 0.183523497750480*I,
	                       -8.238258033274786 + 2.796192505614474*I,
	                       18.785982629476070 - 20.237292093573895*I,
	                      -11.799383335697918 + 46.411650777279597*I,
	                      -11.799383335697890 - 46.411650777279569*I,
	                       18.785982629476067 + 20.237292093573895*I,
	                       -8.238258033274763 - 2.796192505614448*I,
	                        1.319411815998138 - 0.183523497750480*I,
	                       -0.068571505514865 + 0.038419074245888*I,
	                        0.000818433612497 - 0.000581353207069*I};
// Constant used in Caratheodory-Fejer procedure to approximate matrix exponent
const double complex z[] = {-6.998688082445778 - 13.995917029301355*I,
                            -2.235968223749446 - 11.109296400461870*I,
                             0.851707264834878 - 8.503832905826961*I,
                             2.917868800307170 - 6.017345968518187*I,
                             4.206124506834328 - 3.590920783130140*I,
                             4.827493775040721 - 1.193987999180278*I,
                             4.827493775040721 + 1.193987999180278*I,
                             4.206124506834328 + 3.590920783130140*I,
                             2.917868800307170 + 6.017345968518187*I,
                             0.851707264834878 + 8.503832905826961*I,
                            -2.235968223749446 + 11.109296400461870*I,
                            -6.998688082445778 + 13.995917029301355*I};
#else
#pragma GCC error "Incorrect numer of poles for matrix exponentiation"
#endif

SnapperCoreVariables * init_core(int K, double mutation_v, double mutation_u)
{
    SnapperCoreVariables *temp;
    temp = (SnapperCoreVariables * )malloc(sizeof(SnapperCoreVariables));
    //Store number of basis function
	temp-> K = K;

	temp->mutation_v = mutation_v; 
	temp->mutation_u = mutation_u;  
	//Store FFT plans and settings for lfftw3 library 
	temp-> cheby_complex = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2 * (K-1));
	temp-> p_forward = fftw_plan_dft_1d(2*(K-1), temp-> cheby_complex, temp-> cheby_complex,FFTW_FORWARD, FFTW_MEASURE);
	temp-> p_backward = fftw_plan_dft_1d(2*(K-1), temp-> cheby_complex, temp-> cheby_complex,FFTW_BACKWARD, FFTW_MEASURE);
    temp-> root_cheby_complex = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2 * (K-3));
	temp-> p_root_forward = fftw_plan_dft_1d(2*(K-3), temp->root_cheby_complex, temp->root_cheby_complex,FFTW_FORWARD, FFTW_MEASURE);
	temp-> p_root_backward = fftw_plan_dft_1d(2*(K-3), temp->root_cheby_complex , temp->root_cheby_complex,FFTW_BACKWARD, FFTW_MEASURE);
	
	//Store systems for Caratheodory-Fejer procedures 
	temp->A_Q = (double complex ***)malloc(CF_DEG * sizeof(double complex  **));
	for (int p = 0; p < CF_DEG; ++p)
	{
		temp->A_Q[p] = (double complex **)malloc(K * sizeof(double complex *));
		for (int k = 0; k < K; ++k)
			temp->A_Q[p][k] = (double complex *)malloc(K * sizeof(double complex));
	}

	temp->A_Y = (double complex ***)malloc(CF_DEG * sizeof(double complex  **));
	for (int p = 0; p < CF_DEG; ++p)
	{
		temp->A_Y[p] = (double complex **)malloc(K * sizeof(double complex *));
		for (int k = 0; k < K; ++k)
			temp->A_Y[p][k] = (double complex *)malloc(K * sizeof(double complex));
	}

	temp->b_Y = (double *)malloc(K * sizeof(double ));

	return temp;
}

int free_core(SnapperCoreVariables * data)
{
	//BEGIN - Garbage clean-up
	fftw_destroy_plan(data->p_forward);
	fftw_destroy_plan(data->p_backward);
	fftw_free(data->cheby_complex);
	fftw_destroy_plan(data->p_root_forward);
	fftw_destroy_plan(data->p_root_backward);
	fftw_free(data->root_cheby_complex);

	for (int p = 0; p < CF_DEG; ++p)
	{
		for (int r = 0; r < data->K; ++r)
			free(data->A_Q[p][r]);
	}
	for (int p = 0; p < CF_DEG; ++p)
			free(data->A_Q[p]);
	free(data->A_Q);
	
	for (int p = 0; p < CF_DEG; ++p)
	{
		for (int r = 0; r < data->K; ++r)
			free(data->A_Y[p][r]);
	}
	for (int p = 0; p < CF_DEG; ++p)
			free(data->A_Y[p]);
	free(data->A_Y);
	free(data->b_Y);
	// END - Garbage Clean-up
	return 0;
}

// BEGIN - Matrix print functions for debugging purpose
void print_matrix(double ** A, int K)
{
	for (int i = 0; i < K; ++i)
	{
		for (int j = 0; j < K; ++j)
	{
		printf("%.4f,", A[i][j]);
	}
	printf(";\n");
	}	
	printf(";\n");
}

void print_complex_matrix(double complex ** A, int K)
{
	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < K; j++)
	{
		printf("%.4f + %.4fi,", creal(A[i][j]), cimag(A[i][j]));
	}
	printf(";\n");
	}	
	printf(";\n");
}
// END - Print functions for debugging purpose

// Computes  n!
long factorial(int n)
{
  if (n == 0) // Base case
    return 1;
  else
    return (n * factorial(n - 1));
}

//Slow square matrix-matrix multiplicaiton for matrices of size K by K  
//( For test purposes, easy way to set up Matrix operators )
//( Only use this in initilization step so optimal performance not that important ) 
void matrix_multiplication(double ** A, double ** B, double ** C, int K) 
{
    //Initialize C
    for (int i = 0; i < K; ++i) 
    {
        for (int j = 0; j < K; ++j) 
        {
            C[i][j] = 0;
        }
    }
    // Multiplying A and B storing in C.
    for (int i = 0; i < K; ++i) {
        for (int j = 0; j < K; ++j) {
            for (int k = 0; k < K; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}
//Slow matrix-vector multiplicaiton of A (size K by K) and a (K by 1)
//( For test purposes, easy way to set up Matrix operators )
//( Only use this in initilization step so optimal performance not that important )
void matrix_vector_multiplication(double ** A, double *a, double *b, int K )
{
	for ( int i = 0; i < K; i++ )
	{
		b[i] = 0;
		for ( int j = 0; j < K; j++ )
		b[i] += A[i][j] * a[j];
	}
}

//BEGIN - Building blocks for Q matrix  
// ( For test purposes )
// ( Only use this in initilization step so optimal performance not that important )
/*
Matrix M follows from formulas for the product of monomial and a shifted Chebyshev polynomial
See Mason and Handscomb (2002, pg. 101) for derivations 
*/
void matrix_M(double ** M, int K)
{
	M[0][0] = 0.5;
	M[0][1] = 0.25;	
	for( int k = 1; k < K-1; ++k)
	{
		M[k][k] = 0.5;
		M[k][k-1] = 0.25;
		M[k][k+1] = 0.25;
	}	
	M[1][0] = 0.5;
	M[K-1][K-1] = 0.5;
	M[K-1][K-2] = 0.25;

	//print_matrix(M, K);
}
/*
Matrix S follows from formulas for the indefinite integral of shifted Chebyshev polynomial
See Mason and Handscomb (2002, pg. 101) for derivations
*/
void matrix_S(double ** S, int K)
{
	for( int k = 1; k < K - 1; ++k)
	{
		S[k-1][k] = -1 / ( 4 * (double) (k - 1) );
		S[k+1][k] = 1 / ( 4 * (double) (k + 1) );
	}	
	S[0][0] = 0;
	S[0][1] = 0;
	S[1][0] = 0.5;
	S[2][1] = 0.125;
	S[K-2][K-1] = -1 / ( 4 * (double) (K - 2) );
	
	//print_matrix(S, K);
}
/*
Matrix D follows from formulas for the derivative of shifted Chebyshev polynomial
See Mason and Handscomb (2002, pg. 101) for derivations
*/
void matrix_D(double ** D, int K)
{
	for(int j = 0; j < K; ++j )
	{
		if ( j % 2 != 0 )
		{
			D[0][j] = 2.0 * (double) j;
		}
	}	
	for( int j = 0; j < K; ++j)
	{
		for( int i = 0; i < j; ++i)
		{
		if ( i % 2 == 0 && j % 2 == 0 )
		   {
			D[i+1][j] = 4.0 * (double) j;
			}
		
		if ( i % 2 != 0 && j % 2 != 0 )
		   {
			D[i+1][j] = 4.0 * (double) j;
			}
		}
	}

	//print_matrix(D, K);
}
//END - Building blocks for Q matrix

/*
Matrix Q =  \beta_2*D - (beta_1+\beta_2)*S*D + 0.5*S*(1-S)*D^2
follows from the backward diffusion approximation of the Wright-Fisher model
g(x,t) = \beta_2(1-x) - \beta_1 x) \frac{\partial}{\partial x} g(x,t) + \frac{1}{2}x(1-x)\frac{\partial^2}{\partial x^2}
 - See Oksendal (2003, chap.8, pg. 133)
and fromuals for shifted Chebyshev polynomials
 - See Mason and Handscomb (2002, pg. 101) for derivations

!!!Note that solving Q on a species tree we need to divide by \theta (population size)!!! 
*/
void matrix_Q(double ** Q, double beta_1, double beta_2,  int K)
{
	// dynamically create array of pointers of size K
	double **M;
	M = (double **)malloc(K * sizeof(double *));
	
	double **S;
	S = (double **)malloc(K * sizeof(double *));
	
	double **D;
	D = (double **)malloc(K * sizeof(double *));
	
	double **Q_temp;
	Q_temp = (double **)malloc(K * sizeof(double *));
	
	double **Q_temp1;
	Q_temp1 = (double **)malloc(K * sizeof(double *));
	
	for (int r = 0; r < K; ++r)
		M[r] = (double *)malloc(K * sizeof(double));
	
	for (int r = 0; r < K; ++r)
		S[r] = (double *)malloc(K * sizeof(double));
		
	for (int r = 0; r < K; ++r)
		D[r] = (double *)malloc(K * sizeof(double));
		
	for (int r = 0; r < K; ++r)
		Q_temp[r] = (double *)malloc(K * sizeof(double));
	
	for (int r = 0; r < K; ++r)
		Q_temp1[r] = (double *)malloc(K * sizeof(double));
	// Get building blocks for Q
	matrix_M( M, K);
	matrix_S( S, K );
	matrix_D( D, K );	
	//double theta = (beta_1 * beta_2) / (beta_1 + beta_2); 
	// Compute Q = (beta_1) D
	for (int i = 0; i < K; ++i) 
	{
        for (int j = 0; j < K; ++j) 
        {
        	Q[i][j] = beta_1 * D[i][j];
        }
    }		
    
    // Compute Q = (beta_1 + beta_2) D - theta M D
	matrix_multiplication(M, D, Q_temp, K); 
	for (int i = 0; i < K; ++i) 
    {
        for (int j = 0; j < K; ++j) 
        {
            Q[i][j] -= (beta_1 + beta_2) * Q_temp[i][j];
        }
    }
	
	//print_matrix(Q , K);
	// Compute Q = (1/2) D - theta M D + 1/(2)*M*D^2
    matrix_multiplication(D, D, Q_temp, K);
    
    matrix_multiplication(M, Q_temp, Q_temp1, K);
	
	for (int i = 0; i < K; ++i) 
	{
        for (int j = 0; j < K; ++j) 
        {
        	Q[i][j] += 0.25 * Q_temp1[i][j];
        }
    }	
    //print_matrix(Q , K);
	//Compute Q = (theta/2) D - theta M D + 1/(2\theta)*M*(1-M)*D^2
	matrix_multiplication(M, Q_temp1, Q_temp, K);
	
	for (int i = 0; i < K; ++i) 
	{
        for (int j = 0; j < K; ++j) 
        {
        	Q[i][j] -= 0.25 * Q_temp[i][j];
        }
    }	
	
	//matrix_multipliication(S, Q, Q_temp, K);
    //matrix_multipliication(S, Q_temp, Q, K);
    
	// printf("Q matrix \n");
    // print_matrix(Q , K);
	// printf("\n");
	// BEGIN - Garbage Clean-up
	for (int r = 0; r < K; ++r)
		free(M[r]);
	free(M);
	for (int r = 0; r < K; ++r)
		free(S[r]);
	free(S);
	for (int r = 0; r < K; ++r)
		free(D[r]);
	free(D);
	for (int r = 0; r < K; ++r)
		free(Q_temp[r]);
	free(Q_temp);
	for (int r = 0; r < K; ++r)
		free(Q_temp1[r]);
	free(Q_temp1);
	// END - Garbage Clean-up
}

/*
This function sets up systems of the form
{(Q-z_i*I) = v_i, (S2*Q-z_i*S^2) = S^2*v_i}
and stores it in arrays of size (K times K times CF_DEG )
{A_Q, A_Y}
for use in the Caratheody-Fejer procedure to approximate matrix 
exponential exp(Q * T), see (Schmelzer and Trefethen, 2007).
*/
void set_up_systems_for_CF_procedure(double complex *** A_Q, double complex *** A_Y, double * b_Q, double * b_Y, double beta_1, 
										double beta_2,  int K, double delta_t)
{
	// BEGIN - Assign variables	
	double **Q;
	Q = (double **)malloc(K * sizeof(double *));
	for (int r = 0; r < K; ++r)
		Q[r] = (double *)malloc(K * sizeof(double));
	
	double **S;
	S = (double **)malloc(K * sizeof(double *));
	for (int r = 0; r < K; ++r)
		S[r] = (double *)malloc(K * sizeof(double));
	
	double **S2;
	S2 = (double **)malloc(K * sizeof(double *));
	for (int r = 0; r < K; ++r)
		S2[r] = (double *)malloc(K * sizeof(double));
	
	double **A_temp;
	A_temp = (double **)malloc(K * sizeof(double *));
	for (int r = 0; r < K; ++r)
		A_temp[r] = (double *)malloc(K * sizeof(double));
	// END - Assign variables
 
	matrix_Q( Q, beta_1 , beta_2, K);
	matrix_S( S, K);	
	for (int p = 0; p < CF_DEG; ++p)
	{
		//A_Q[][][p] = Q - z[p]*I
		for (int i = 0; i < K; ++i) 
		{
		    for (int j = 0; j < K; ++j) 
		    {
		    	A_Q[p][i][j] = Q[i][j];
		    	if ( i == j)
		    	{
		    		A_Q[p][i][j] -= (z[p] / delta_t);
		    	}	    	
		    }
		}
		// printf("A_Q %d", p);
		// printf("\n");
		// print_complex_matrix(A_Q[p], K);
				
		//A_Y[][][p] = S^2Q - z[p]S^2
		matrix_multiplication(S,S,S2, K);
		matrix_multiplication(S2, Q, A_temp, K);	
		for (int i = 0; i < K; ++i) 
		{
		    for (int j = 0; j < K; ++j) 
		    {
		    	A_Y[p][i][j] = A_temp[i][j] - (z[p] / delta_t) * S2[i][j]  ;
		    }
		}
		// printf("A_Y %d", p);
		// printf("\n");
		// print_complex_matrix(A_Y[p], K);
	}
	
	//b_Y = S2*b_Q

	for (int i = 0; i < K; ++i) 
		    {
		    	b_Q[i] = b_Q[i] / delta_t  ;
		    }

	matrix_vector_multiplication( S2, b_Q, b_Y, K);
	
	//printf("b_Y");
	//printf("\n");
	//for (int r = 0; r < K; ++r)
	//	printf("%.6f,",b_Y[r]);
	
	// BEGIN - Garbage Clean-up
	for (int r = 0; r < K; ++r)
		free(Q[r]);
	free(Q);
	for (int r = 0; r < K; ++r)
		free(S[r]);
	free(S);
	for (int r = 0; r < K; ++r)
		free(A_temp[r]);
	free(A_temp);
	for (int r = 0; r < K; ++r)
		free(S2[r]);
	free(S2);
	// END - Garbage Clean-up
}

int set_matrix_Q(double * lambda_0, double theta, double t, SnapperCoreVariables * tree)
{
    double beta_1 = theta * ( tree->mutation_u / (tree->mutation_u + tree->mutation_v) );
	double beta_2 = theta * ( tree->mutation_u / (tree->mutation_u + tree->mutation_v) );
	set_up_systems_for_CF_procedure(tree->A_Q, tree->A_Y, lambda_0, tree->b_Y, beta_1, beta_2, tree->K, t);
    return 0;
}

/*
Solves linear systems in Caratheodory-Fejer procedure 
when \beta_1 = beta_2
using  A_Q x = b_Q and A_Y x = b_Y
where b_Q = \lambda(0) and b_Y = S^2 \lambda(0)
(Note that these two systems are equivalent) 
*/
void sparse_solver_CF(double complex ** A_Q, double complex ** A_Y, double * b_Q, double * b_Y, double complex * x,   int K)
{
// Gaussian backward substitution (for banded matrix with every second matrix  
// element in band equal to zero and bandwidth = 5) 
// to solve x using the 
// two systems A_Q x = b_Q and A_Y x = b_Y
	x[K-1] = b_Q[K-1] / A_Q[K-1][K-1];
    x[K-2] = b_Q[K-2] / A_Q[K-2][K-2];
    x[K-3] = (b_Y[K-1] - A_Y[K-1][K-1] * x[K-1]) / A_Y[K-1][K-3]; 
    x[K-4] = (b_Y[K-2] - A_Y[K-2][K-2] * x[K-2]) / A_Y[K-2][K-4];
    for (int i = K - 5; i > -1 ; --i)
		x[i] = (b_Y[i+2] - A_Y[i+2][i+2] * x[i+2] - A_Y[i+2][i+4] * x[i+4]) / A_Y[i+2][i];
        
}
/*
Solves linear system in Caratheodory-Fejer procedure 
when \beta_1 != beta_2
based on  A_Q x = b_Q and A_Y x = b_Y 
where b_Q = \lambda(0) and b_Y = S^2 \lambda(0)
(Note that these two systems are equivalent) 
*/
void solver_CF(double complex ** A_Q, double complex ** A_Y, double * b_Q, double * b_Y, double complex * x,   int K)
{
// Gaussian backward substitution (for banded matrix w. bandwidth = 5) to solve x but using the 
// two systems A_Q x = b_Q and A_Y x = b_Y.
	x[K-1] = b_Q[K-1] / A_Q[K-1][K-1];
    x[K-2] = (b_Q[K-2] - A_Q[K - 2][K-1] * x[K-1]) / A_Q[K-2][K-2];
    x[K-3] = (b_Y[K-1] - A_Y[K-1][K-1] * x[K-1] 
    		- A_Y[K-1][K-2] * x[K-2]) / A_Y[K-1][K-3];
    x[K - 4] = (b_Y[K-2] - A_Y[K-2][K-1] * x[K-1] - A_Y[K-2][K-2] * x[K-2] 
    		- A_Y[K-2][K-3] * x[K-3]) / A_Y[K-2][K-4];
    for (int i = K - 5; i > -1; --i)
	{	
		x[i] = b_Y[i+2]; 
		for (int j = 1; j < 5; ++j)
		{
			x[i] -= A_Y[i+2][i+j] * x[i+j];
		}
    x[i] /= A_Y[i+2][i];
    }
   
}
/*
Solve linear system of Chebyshev coefficients at the root of tree,
for g(x) given Chebyshev coefficients of f(x) 
such that f(x) = x(1-x)g(x) holds.
*/
void solver_root(double *coeff, int K)
{
	coeff[K-1] = -16.0 * coeff[K-1];
	coeff[K-2] = -16.0 * coeff[K-2];	
	coeff[K-3] = -16.0 * (-0.125 * coeff[K-1] + coeff[K-3]);
	coeff[K-4] = -16.0 * (-0.125 * coeff[K-2] + coeff[K-4]);
	for (int  k = K - 5; k > 2; --k) 
	{
		coeff[k] = -16.0 * (-0.125 * coeff[k+2] + 0.0625*coeff[k+4] + coeff[k]); 
	}
	coeff[2] = -8.0 * (-0.125 * coeff[4] + 0.0625*coeff[6] + coeff[2]);
/*g(x) is a K-2 Chebyshev basis approximation*/
	coeff[1] = 0;
	coeff[0] = 0;	
}
/*
Caratheodory-Fejer procedure to approximate matrix 
exponential exp(Q * t), see (Schmelzer and Trefethen, 2007).
*/
int matrix_exponent(double complex *** A_Q, double complex *** A_Y, double * lambda_0, double * b_Y, double * lambda_t, double mutation_u, double mutation_v, double t, double step_size, int K)
{
    double complex * x;
	x = (double complex *)malloc(K * sizeof(double complex)); 
	//Note that we always assume t = stepsize

	//printf("%g, %g,", t, step_size);
	for( int k = 0; k < K; ++k)
			{ 
				lambda_t[k] =  0; 
			}

	for( double step = 0; step < t; step+= step_size )
	{
	//printf("%g;", step);
		for( int p = 0; p < CF_DEG; ++p)
		{
	//			printf("%d,", i);
    //Banded matrices,A_Q and A_Y, have zero entries when \beta_1 = \beta_2
            if((mutation_u - mutation_v) > 0)
			{			
			solver_CF(A_Q[p], A_Y[p], lambda_0, b_Y, x,  K);
			}
			else
			{			
			sparse_solver_CF(A_Q[p], A_Y[p], lambda_0, b_Y, x,  K);		
			}
			
			for( int k = 0; k < K; ++k)
			{ 
	//lambda_t is the Chebyshev coefficients at top of a branch with length t.  
				// printf("%d,", k);
				// printf("%.4f + %.4fi,", creal(x[k]), cimag(x[k]));
    			// printf("\n");
				lambda_t[k] += creal(x[k]*c[p]); 
			}
		}

	}	
    //BEGIN - Garbage clean-up
    free(x);
    //END - Garbage clean-up
	return 0;
}
//Userfriendly wrapper for the matrix exponent function
double * calculate_lambda_t(double * lambda_0,  double t,  SnapperCoreVariables * tree)
{	
	double * lambda_t;
    lambda_t = (double *)malloc(tree->K * sizeof(double ));

	matrix_exponent(tree->A_Q, tree->A_Y, lambda_0, tree->b_Y, lambda_t, tree->mutation_u, tree->mutation_v, t, t, tree->K);

	return lambda_t;
}
/*
Integrate function on domain [0,1] given function K Chebyshev coefficients
using Clenshaw-Curtis type quadrature (see Trefethen, 2013, pg 192) 
*/
double Clenshaw_Curtis_integration(double *coeff, int K)
{
	double integral = 0;
	for (int  k = 0; k < K; ++k) 
	{
		if (k % 2 == 0)
        	integral += coeff[k] / (1 - pow(k,2));
	}
	return integral;
}
/*
Computes 
\int_0^1 f(x)\pi(\beta_1,\beta2)dx

where f(x) = x(1-x)g(x)

using Clenshaw-Curtis quadrature.

This gives the likelihood of a tree at a particular site given the combined partial likelihood coefficients at the bottom of 
the root branch.
*/
double integrate_at_root(double * partial_likelihood_bottom, double theta, SnapperCoreVariables * tree)
{	
	double x[tree->K-2]; 
	double root_integral = 0;
	
	double f_0 = partial_likelihood_bottom[0];
	double f_1 = partial_likelihood_bottom[tree->K-1];
	double beta_1 = theta * (tree->mutation_u / ( tree->mutation_u + tree->mutation_v));
	double beta_2 = theta * (tree->mutation_u / ( tree->mutation_u + tree->mutation_v));

	transform_to_chebyshev_coef(partial_likelihood_bottom, tree);

	solver_root(partial_likelihood_bottom, tree->K); 
	//printf("Backwards substitution \n");
	//for (int i=0; i < K; ++i){
	//	printf("%g,", cheby_values[i]);
	//}
	//printf("\n\n");
	
	transform_to_chebyshev_values_root(&partial_likelihood_bottom[2], tree->root_cheby_complex, tree->p_root_forward, tree->K-2);
	generate_clobatto_grid(x,tree->K-2);
	for (int k = 0; k < tree->K-2; ++k)
	{
		partial_likelihood_bottom[k] = pow(x[k], beta_1) * pow((1 - x[k]), beta_2) * partial_likelihood_bottom[2+k];
	}
	transform_to_chebyshev_coef_root(partial_likelihood_bottom, tree->root_cheby_complex, tree->p_root_backward, tree->K-2);;
	
	root_integral += exp(lgamma(beta_1 + beta_2) 
                    - lgamma(beta_1) - lgamma(beta_2))
		*Clenshaw_Curtis_integration(partial_likelihood_bottom, tree->K-2);
	
	root_integral += f_0 + beta_1 / (beta_1 + beta_2) * (f_1 - f_0); 
	
	return root_integral;
}
/*
Evaluate partial likelihood at bottom of leave branch on grid x (of size K) 
given input data n (sample size), r (number of red alleles in sample)
*/
double * fit_leaf_likelihood_function_values(int n, int r, SnapperCoreVariables * tree) 
{
	double * x, * function_values;

	x = (double *)malloc(tree->K * sizeof(double ));
	function_values = (double *)malloc(tree->K * sizeof(double ));

	generate_clobatto_grid(x, tree->K);

	for (int k = 0; k < tree->K; ++k) 
	{
		function_values[k] = ( factorial(n) / (factorial(r) * factorial(n - r)) )
						   * pow(x[k],r) * pow(1 - x[k],(n-r));
	}
	free(x);
	return function_values;
}
/*
Generate Chebyshev-Lobatto grid with K points (K equispaced points on unit semi-circle
mapped to [0,1], see Trefethen 2013) 
*/
void generate_clobatto_grid(double *x, int K) 
{   
	for (int  k = 0; k < K; ++k) 
	{
		x[k] = 0.5 - cos( -k / (K-1.0) * M_PI) / 2.0;
	}
}
/*
Use the backward FFT algorithm to get K Chebyshev coefficients given partial likelihood evaluated on Chebyshev-Lobatto grid
with K points    
*/
void transform_to_chebyshev_coef(double *cheby_values, SnapperCoreVariables * tree) 
{
	for (int k = 0; k < tree->K; ++k)
	{
		tree->cheby_complex[k] = cheby_values[tree->K-1-k];
	}
	for (int k = 1; k < tree->K-1; ++k)
	{
		tree->cheby_complex[tree->K-1 + k] = cheby_values[k];
	}
	fftw_execute(tree->p_backward);
	for (int k = 0; k < tree->K; ++k)
	{
		cheby_values[k] = creal(tree->cheby_complex[k]) / (double) (tree->K-1);
	}
	cheby_values[0] = cheby_values[0] / 2.0;
}
/*
Use the forward FFT algorithm to get partial likelihood evaluated on Chebyshev-Lobatto grid
with K points given K Chebyshev coefficients     
*/
void transform_to_chebyshev_values(double *cheby_coeff, SnapperCoreVariables * tree)
{	
	
	for (int k = 0; k < tree->K; ++k)
	{
		cheby_coeff[k] = cheby_coeff[k] * (double) (tree->K-1);
	}
		cheby_coeff[0] = cheby_coeff[0] * 2.0;
	
	for (int k = 0; k < tree->K; ++k)
	{
		tree->cheby_complex[k] = cheby_coeff[k];
	}
	for (int k = 1; k < tree->K-1; ++k)
	{
		tree->cheby_complex[tree->K-1 +k] = cheby_coeff[tree->K-1-k];
	}
	
	fftw_execute(tree->p_forward);
	for (int k = 0; k < tree->K; ++k)
	{
		cheby_coeff[k] = creal(tree->cheby_complex[tree->K-1 - k]) / (2 * (tree->K-1));
	}
}
// Use K-2 number of basis functions at root therefore we use different FFT plans   
void transform_to_chebyshev_coef_root(double *cheby_values,  fftw_complex *cheby_complex, fftw_plan plan, int K)
{
	for (int k=0; k < K; ++k)
	{
		cheby_complex[k] = cheby_values[K-1-k];
	}
	for (int k=1; k < K-1; ++k)
	{
		cheby_complex[K-1 +k] = cheby_values[k];
	}
	fftw_execute(plan);
	for (int k=0; k < K; ++k)
	{
		cheby_values[k] = creal(cheby_complex[k]) / (double) (K-1);
	}
	cheby_values[0] = cheby_values[0] / 2.0;
}
// Use K-2 number of basis functions at root therefore we use different FFT plans 
void transform_to_chebyshev_values_root(double *cheby_values, fftw_complex *cheby_complex, fftw_plan plan, int K)
{	
	
	for (int k=0; k < K; ++k)
	{
		cheby_values[k] = cheby_values[k] * (double) (K-1);
	}
		cheby_values[0] = cheby_values[0] * 2.0;
	
	for (int k=0; k < K; ++k)
	{
		cheby_complex[k] = cheby_values[k];
	}
	for (int k=1; k < K-1; ++k)
	{
		cheby_complex[K-1 +k] = cheby_values[K-1-k];
	}
	
	fftw_execute(plan);
	for (int k=0; k < K; ++k)
	{
		cheby_values[k] = creal(cheby_complex[K-1-k])/(2*(K-1));
	}
}

