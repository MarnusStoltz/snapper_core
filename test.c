#include "snapper_core.h"

int main() 
{
    SnapperCoreVariables *tree;
    // SETTINGS
    int K = 7; //Number of Chebyshev basis funcitons

    // INPUT DATA FOR TREE (with two taxa) 
    int n_1 = 10; //Total number of samples at site for species 1
    int r_1 = 3; //Number of red alleles observed at site for species 1
    double theta_1 = 2.0; //Population size on leaf branch 1
    double t_1 = 0.1; //Branch length

    int n_2 = 10; //total number of samples at site for species 2
    int r_2 = 3; //number of red alleles observed at site for species 2
    double theta_2 = 2.0; //Population size on leaf branch 2
    double t_2 = 0.1; //Branch length
    
    double theta_r = 0.1; //Population size on root branch

    // VARIABLES
    double * partial_likelihood_bottom; // Store partial likelihood along a branch
    double * partial_likelihood_top; // Store partial likelihood along a branch
    double likelihood;

    //Initialize core variables and FFT plans 
    tree = init_core(K,1.0,1.0);

    //LEAF BRANCH 1
    partial_likelihood_bottom = fit_leaf_likelihood_function_values(n_1, r_1, tree);

    printf("Function values on Cheby-Lobatto grid (bottom of branch 1) \n");
	for (int i = 0; i < K; ++i){
		printf("%g,", partial_likelihood_bottom[i]);
	}
    printf("\n\n");

    transform_to_chebyshev_coef(partial_likelihood_bottom, tree);
    
    printf("Chebyshev coefficients (bottom of branch 1) \n");
	for (int i = 0; i < K; ++i){
		printf("%g,", partial_likelihood_bottom[i]);
	}
    printf("\n\n");

    //Setup backward PDE in terms of Chebyshev basis functions
    set_matrix_Q(partial_likelihood_bottom, theta_1, t_1 , tree);
    
    //Get Chebyshev coeffients of partial likielihood at top of branch by solving PDE along the branch for length t_0, 
    partial_likelihood_top = calculate_lambda_t(partial_likelihood_bottom, t_1, tree);

    printf("Chebyshev coefficient (top of branch 1) \n");
    
	for (int i = 0; i < K; ++i){
		printf("%g,", partial_likelihood_top[i]);
	}
    printf("\n\n");

    //Get function values of partial likelihood at top of branch at Chehbyshev Lobatto points
    transform_to_chebyshev_values(partial_likelihood_top, tree); 

    printf("Function values on Cheby-Lobatto grid (top of branch 1) \n");
	for (int i = 0; i < K; ++i){
		printf("%g,", partial_likelihood_top[i]);
	}
    printf("\n\n");
    //LEAF BRANCH 2
    //Note that input settings are identical therefore don't compute anything; just use partial likelihood 
    //already calculated for LEAF BRANCH 1 

    //ROOT BRANCH 
    //Get chebyshev function values for partial likelihood at bottom of root branch
	for (int i = 0; i < K; ++i){
		partial_likelihood_bottom[i] = partial_likelihood_top[i]*partial_likelihood_top[i];
	}
    //Compute the likelihood of a tree by integrating over the partial likelihood at bottom of 
    //root branch multiplied with with stationary beta distribution 
    likelihood = integrate_at_root(partial_likelihood_bottom, theta_r, tree);
    printf("Site likelihood of tree\n");
    printf("%g,", likelihood);
    printf("\n\n");

    //BEGIN - Garbage clean-up
    free_core(tree);
    //END - Garbage clean-up
	return 0;
}
