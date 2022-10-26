#### SIGMA ####
posterior_sigma <- function(sigma, alpha_sigma, beta_sigma,
                            theta, alpha_theta, beta_theta,
                            rho){
  
  
  k <- length(rho)
  
  posterior <- (alpha_sigma - 1) * log(sigma) + (beta_sigma - 1) * log(1 - sigma)
  
  posterior <- posterior + (alpha_theta - 1) * log(theta + sigma) + sigma * (- beta_theta) 
  
  for(i in 1:(k-1)){
    
    n_clust <- rho[i]
    posterior <- posterior + log(theta + i * sigma) + log_pochhammer(1 - sigma, n_clust - 1)
    
    
  }
  
  
  posterior <- posterior + log_pochhammer(1 - sigma, rho[k])
  
  return(posterior)
  
  
  
  
}



