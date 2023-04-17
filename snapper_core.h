//#ifndef SNAPPER_CORE_H
//#define SNAPPER_CORE_H

#include <stdlib.h>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <complex.h>
#include <cblas.h>
#include <fftw3.h>

//Number of poles used in Caratheodory Fejer procedure
#define CF_DEG 12

//STRUCTURES
//Store settings that will remain fixed for a tree
struct SnapperCore{
    
    //Number of Chebyshev basis funciton
    int K;
    //Forward and backward mutation ratios
    double mutation_u;
    double mutation_v;
    //Optimised FFT algorithm execution plan for specific machine
    fftw_complex * cheby_complex;
	fftw_plan p_forward, p_backward;
    fftw_complex *root_cheby_complex;
	fftw_plan p_root_forward, p_root_backward;
	//Pointers for systems used in Caratheodory-process to approximate matrix exponential
    double complex *** A_Q;	
	double complex *** A_Y;
    double * b_Y;


};
typedef struct SnapperCore SnapperCoreVariables;

//LIST OF FUNCTIONS
long factorial(int );
void matrix_multiplication(double **, double **, double **, int );
void matrix_vector_multiplication(double **, double *, double *, int );
void print_matrix(double **, int);
void print_complex_matrix(double complex **, int);
SnapperCoreVariables * init_core(int , double, double);
int free_core(SnapperCoreVariables *);
double * fit_leaf_likelihood_function_values(int , int , SnapperCoreVariables *);
int set_matrix_Q(double *, double , double , SnapperCoreVariables * );
void matrix_D(double **, int );
void matrix_M(double **, int );
void matrix_S(double **, int );
void matrix_Q(double **, double, double, int );
void set_up_systems_for_CF_procedure(double complex ***, double complex ***, double *, double *, 
                                        double, double,  int , double);
void sparse_solver_CF(double complex **, double complex **, double *, double *, double complex *,   int );
void solver_CF(double complex **, double complex **, double *, double *, double complex *,   int );
void solver_root(double *, int );
int matrix_exponent(double complex ***, double complex ***, double *, double *, double *, double, double, double, double, int );
double * calculate_lambda_t(double *,  double,  SnapperCoreVariables * );
double Clenshaw_Curtis_integration(double *, int );
double integrate_at_root(double * , double, SnapperCoreVariables * );
void generate_clobatto_grid(double *, int ); 
void transform_to_chebyshev_coef(double *, SnapperCoreVariables *); 
void transform_to_chebyshev_values(double * , SnapperCoreVariables * );
void transform_to_chebyshev_coef_root(double *,  fftw_complex *, fftw_plan , int);
void transform_to_chebyshev_values_root(double *, fftw_complex *, fftw_plan, int);



