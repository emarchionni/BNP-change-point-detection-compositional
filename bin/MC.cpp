#include <Rcpp.h>

// [[Rcpp::export]]
double MC_log_integrated_MARGINAL_likelihood_cluster_single_cpp(Rcpp::NumericVector y, int trunc, int iter_omega, double alpha_omega, double beta_omega){
	
	double value = .0;
	int d = y.length();
	Rcpp::NumericVector omega_sampled(d);
	
	
	for(int iter = 0; iter < iter_omega; ++iter){
		
		omega_sampled = Rcpp::rgamma(d, alpha_omega, 1/beta_omega );
		value += std::exp(log_ddirichlet_cpp(y, omega_sampled);
		
		
	}
	
	
	return std::log(value/iter_omega);
	
}

// [[Rcpp::export]]
double MC_log_integrated_MARGINAL_likelihood_cluster_multiple_cpp(Rcpp::NumericMatrix y, int trunc, int iter_omega, double alpha_omega, double beta_omega){
	
	double value = .0;
	int d = y.ncol();
	Rcpp::NumericVector omega_sampled(d);
	
	for(int iter = 0; iter < iter_omega; ++iter){
		
		omega_sampled = Rcpp::rgamma(d, alpha_omega, 1/beta_omega );
		value += log_integrated_likelihood_cluster_multiple_cpp(y, omega_sampled, trunc);
		
	}
	
	return std::log(value/iter_omega);
		
}