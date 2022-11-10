#include <Rcpp.h>


// VARIOUS SUPPORT FUNCTIONS


// [[Rcpp::export]]
double pochhammer_cpp(double x, unsigned factor){
  
	double value = 1.;
	
	if(factor == 0)
		return value;
	
	for(unsigned i = 0; i < factor; ++i)
		value *= (x + i);
	
	return value;
}


// [[Rcpp::export]]
double factorial_cpp(double n){
  
  
  double fact_n = 1.;
  
  if(n == 0)
    return fact_n;
  
  for(int i = 1; i <= n; ++i)
    fact_n *= i;
  
  return fact_n;
}


// [[Rcpp::export]]
double prod_factorial_vector_cpp(Rcpp::NumericVector y){
  
  double value = 1.;
  
  for(int i = 0; i < y.length(); ++i)
    value *= factorial_cpp(y[i]);
  
  return value;
  
}


// [[Rcpp::export]]
double sum_cpp(Rcpp::NumericVector y){
	
	double value = 0;
	
	for(int i = 0; i < y.length(); ++i)
		value += y[i];
	
	return value;
}


//// [[Rcpp::export]]	
//double prod_cpp(Rcpp::NumericVector y){
//	
//	double value = 1;
//	
//	for(int i = 0; i < y.length(); ++i)
//		value *= y[i];
//	
//	return value;
//}	

//// [[Rcpp::export]]	
//Rcpp::NumericVector sum_vectors_cpp(Rcpp::NumericVector a, Rcpp::NumericVector b){
//  
//  Rcpp::NumericVector sum(a.length());
//  
//  for(int i = 0; i < a.length(); ++i)
//    sum[i] = a[i] + b[i];
//  
//  return sum;
//  
//}


// [[Rcpp::export]]	
double log_ddirichlet_cpp(Rcpp::NumericVector y, Rcpp::NumericVector omega){
	
	double value = 0.;
	
	for(int i = 0; i < y.length(); ++i){
		value += ( (omega[i] - 1) * std::log(y[i]) - std::lgamma(omega[i]));
	}
	
	value += std::lgamma(Rcpp::sum(omega));
	
	return value;
	
}


// GET INTEGER WEAK COMPOSITIONS


// [[Rcpp::export]]
Rcpp::NumericVector get_first_weak_composition_cpp(int m, int d){

	Rcpp::NumericVector composition(d);
	
	composition[d - 1] = m;
	
	return composition;

}




// [[Rcpp::export]]
bool exist_next_weak_composition_cpp(int m, int d, Rcpp::NumericVector composition){
	
	if(composition[0] == m)
		return false;
	
	return true;
	
	
}


// [[Rcpp::export]]
Rcpp::NumericVector get_next_weak_composition_cpp(int m, int d, Rcpp::NumericVector composition){
	
	int last = d - 1;
	
	while(composition[last] == 0)
		last--;
	
	int z = composition[last];
	
	Rcpp::NumericVector next_composition(composition);

	
	    
    next_composition[last - 1] += 1;
    next_composition[last] = 0;
    next_composition[d - 1] = z - 1;
	
	return next_composition;
}






// [[Rcpp::export]]
Rcpp::NumericMatrix get_all_weak_composition_cpp(int m, int d){
	
	int num_compositions = R::choose(m+d-1, m);
	
	Rcpp::NumericMatrix all_composition(num_compositions, d);
	
	int index_composition = 0;
	
	Rcpp::NumericVector composition = get_first_weak_composition_cpp(m, d);
	
	Rcpp::NumericMatrix::Row row = all_composition.row(index_composition);
	row = composition;
	
	while(exist_next_weak_composition_cpp(m, d, composition)){
		
		++index_composition;
	  
		composition = get_next_weak_composition_cpp(m, d, composition);
		
		all_composition.row(index_composition) = composition;
		
		
	}
	
	return all_composition;
	
}


// TRANSITION DENSITIES

// [[Rcpp::export]]
double zeta_m_cpp(Rcpp::NumericVector y_0, Rcpp::NumericVector y, Rcpp::NumericVector omega, int m){
	
	if(m == 0)
		return 1.;
	
	int d = y.length();
	
	double value = .0;
	
	Rcpp::NumericMatrix compositions = get_all_weak_composition_cpp(m, d);
	
	int num_compositions = compositions.nrow();
	
	for(int i = 0; i < num_compositions; ++i){
		
		Rcpp::NumericVector curr_composition(d);
	  
	  for(int j = 0; j < d; ++j)
	    curr_composition[j] = compositions(i, j);
	  
	  Rcpp::NumericVector sum_omega = omega + curr_composition;
		
		value += ( factorial_cpp(m) * pow(Rcpp::sum(y_0), curr_composition[d - 1]) * std::exp(log_ddirichlet_cpp(y, sum_omega)) / prod_factorial_vector_cpp(curr_composition));
		
	}
	
	return value;
	
}



// [[Rcpp::export]]
double Q_n_poly_cpp(Rcpp::NumericVector y_0, Rcpp::NumericVector y, Rcpp::NumericVector omega, int n){
	
	if(n == 0)
		return 1.;
	
	double omega_norm = Rcpp::sum(omega);
	double value = .0;
	double p = .0;
	double z_m = .0;
	
	for(int m = 0; m <= n; ++m){
		p = pochhammer_cpp(omega_norm + m, n - 1);
		z_m = zeta_m_cpp(y_0, y , omega, m);

		
		value += ( pow(-1, n - m) * R::choose(n, m) * p * z_m );
	}
	

	
	return (omega_norm + 2 * n - 1) * value / factorial_cpp(n);
	
}


// [[Rcpp::export]]
double log_transition_densities(Rcpp::NumericVector y_0, Rcpp::NumericVector y, Rcpp::NumericVector omega, int trunc){
  
  double value = std::exp(log_ddirichlet_cpp(y, omega));
  double Qn = 0;
  double lambda_n = 0;
  double omega_norm = Rcpp::sum(omega);
  
  for(int n = 0; n <= trunc; ++n){
	  
    Qn = Q_n_poly_cpp(y_0, y, omega, n);
    lambda_n = (.5) * n * (n - 1 + omega_norm);
    value += (std::exp(-lambda_n) * Qn);
	
  }
  
  if(value <= 0){
    Rcpp::warning("transition negative or null");
    return 0.;
  }
	  

  return std::log(value);
  
}



// LOG INTEGRATED LIKELIHOOD CLUSTER

// [[Rcpp::export]]
double log_integrated_likelihood_cluster_multiple_cpp(Rcpp::NumericMatrix y, Rcpp::NumericVector omega, int trunc){
	// y: observation in the cluster / dim: time x components
	// omega: omega of cluster / dim: d
	
	int d = omega.length();
	int n_clust = y.nrow();
	Rcpp::NumericVector y_curr(d);
	
	for(int j = 0; j < d; ++j)
		y_curr[j] = y(0, j);
		
	
	double value = log_ddirichlet_cpp(y_curr, omega);
	
	
	Rcpp::NumericVector y_previous(d);
	
	for(int i = 1; i < n_clust; ++i){
		
		y_previous = y_curr;
		
		for(int j = 0; j < d; ++j)
			y_curr[j] = y(i, j);
		
		
		value += log_transition_densities(y_previous, y_curr, omega, trunc);
		
		
	}
	
	return value;
	
}


// MH CLUSTER OMEGA

// [[Rcpp::export]]
Rcpp::NumericVector MH_omega_cluster_multi(int iter_omega, Rcpp::NumericMatrix y, int d, double sigma_0, double alpha_omega, double beta_omega, int trunc){
	
	//Rcpp::NumericMatrix sigma = Rcpp::NumericMatrix::diag(d, sigma_0);
	
	Rcpp::NumericVector omega(d);
	Rcpp::NumericVector proposed_omega(d);
	Rcpp::NumericVector exp_omega(Rcpp::exp(omega));
	Rcpp::NumericVector exp_proposed_omega(Rcpp::exp(proposed_omega));
	
	
	double numerator = .0;
	double denominator = .0;
	double ratio = .0;


	
	for(int iter = 0; iter < iter_omega; ++iter){
		
		
		
		for(int i = 0; i < d; ++i)
			proposed_omega[i] = R::rnorm(omega[i], sigma_0);
		
		exp_proposed_omega = Rcpp::exp(proposed_omega);
		
		numerator = log_integrated_likelihood_cluster_multiple_cpp(y, exp_proposed_omega, trunc);
		denominator = log_integrated_likelihood_cluster_multiple_cpp(y, exp_omega, trunc);
		
		for(int j = 0; j < d; ++j){
			
			// Rcpp::dgamma is shape-scale, our convention is shape-rate
			numerator += R::dgamma(exp_proposed_omega[j], alpha_omega, 1 / beta_omega, 1);
			denominator += R::dgamma(exp_omega[j], alpha_omega, 1 / beta_omega, 1);
			
		}
		
		numerator += Rcpp::sum(proposed_omega);
		denominator += Rcpp::sum(omega);
		
		ratio = numerator - denominator;
		
		if(std::log(R::runif(0.0,1.0)) <= std::min(0.0, ratio)){
			
			omega = proposed_omega;
			exp_omega = exp_proposed_omega;		
			
		}
		
	}
	
	return(exp_omega);
	
}


// [[Rcpp::export]]
Rcpp::NumericVector MH_omega_cluster_single(int iter_omega, Rcpp::NumericVector y, int d, double sigma_0, double alpha_omega, double beta_omega, int trunc){
	
	//Rcpp::NumericMatrix sigma = Rcpp::NumericMatrix::diag(d, sigma_0);
	
	Rcpp::NumericVector omega(d);
	Rcpp::NumericVector proposed_omega(d);
	Rcpp::NumericVector exp_omega(Rcpp::exp(omega));
	Rcpp::NumericVector exp_proposed_omega(Rcpp::exp(proposed_omega));
	
	
	double numerator = .0;
	double denominator = .0;
	double ratio = .0;


	
	for(int iter = 0; iter < iter_omega; ++iter){
		
		
		for(int i = 0; i < d; ++i)
			proposed_omega[i] = R::rnorm(omega[i], sigma_0);
		
		exp_proposed_omega = Rcpp::exp(proposed_omega);
		
		numerator = log_ddirichlet_cpp(y, exp_proposed_omega);
		denominator = log_ddirichlet_cpp(y, exp_omega);
		
		for(int j = 0; j < d; ++j){
			
			// Rcpp::dgamma is shape-scale, our convention is shape-rate
			numerator += R::dgamma(exp_proposed_omega[j], alpha_omega, 1 / beta_omega, 1);
			denominator += R::dgamma(exp_omega[j], alpha_omega, 1 / beta_omega, 1);
			
		}
		
		numerator += Rcpp::sum(proposed_omega);
		denominator += Rcpp::sum(omega);
		
		ratio = numerator - denominator;
		
		if(std::log(R::runif(0.0,1.0)) <= std::min(0.0, ratio)){
			
			omega = proposed_omega;
			exp_omega = exp_proposed_omega;		
			
		}
		
	}
	
	return(exp_omega);
	
}


// MC INTEGRATION LIKELIHOOD

// [[Rcpp::export]]
double MC_log_integrated_MARGINAL_likelihood_cluster_single_cpp(Rcpp::NumericVector y, int trunc, int iter_omega, double alpha_omega, double beta_omega){
	
	double value = .0;
	int d = y.length();
	Rcpp::NumericVector omega_sampled(d);
	
	
	for(int iter = 0; iter < iter_omega; ++iter){
		
		omega_sampled = Rcpp::rgamma(d, alpha_omega, 1/beta_omega );
		value += std::exp(log_ddirichlet_cpp(y, omega_sampled));
		
		
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
		value += std::exp(log_integrated_likelihood_cluster_multiple_cpp(y, omega_sampled, trunc));
		
	}
	
	return std::log(value/iter_omega);
		
}
