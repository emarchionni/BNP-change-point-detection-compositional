#include <Rcpp.h>


// [[Rcpp::export]]
double factorial_cpp(double n){
  
  
  double fact_n = 1;
  
  if(n == 0)
    return fact_n;
  
  for(int i = 1; i <= n; ++i)
    fact_n *= i;
  
  return fact_n;
}


// [[Rcpp::export]]
double prod_factorial_vector_cpp(Rcpp::NumericVector y){
  
  double value = 1;
  
  for(int i = 0; i < y.length(); ++i)
    value *= factorial_cpp(y[i]);
  
  return value;
  
}


// [[Rcpp::export]]
double prod_factorial_vector_cpp_2(Rcpp::NumericVector y){
  
  Rcpp::NumericVector y_fact = Rcpp::factorial(y);
  
  double value = 1.;
  
  for(int i = 0; i < y.length(); ++i)
    value *= y_fact[i];
  
  return value;
  
}
