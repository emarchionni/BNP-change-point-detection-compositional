#### ALPHA SIGMA ####


MC_log_alpha_sigma <- function(sigma_old, sigma_proposed,
                               alpha_sigma, beta_sigma,
                               theta, alpha_theta, beta_theta,
                               rho){
  
  
  log_ratio <- posterior_sigma(sigma_proposed, alpha_sigma, beta_sigma,
                               theta, alpha_theta, beta_theta,
                               rho)
  
  log_ratio <- log_ratio - posterior_sigma(sigma, alpha_sigma, beta_sigma,
                                           theta, alpha_theta, beta_theta,
                                           rho)
  
  return(log_ratio)
  
  
}


#### POSTERIOR SIGMA ####

posterior_sigma <- function(sigma, alpha_sigma, beta_sigma,
                            theta, alpha_theta, beta_theta,
                            rho){
  
  n <- sum(rho)
  k <- length(rho)
  
  posterior <- (alpha_sigma - 1) * log(sigma) + (beta_sigma - 1) * log(1 - sigma)
  
  posterior <- posterior + (alpha_theta - 1) * log(theta + sigma) + sigma * (- beta_theta) 
  
  if(k == 1){
    posterior <- posterior + log(theta + sigma) + log_pochhammer(1 - sigma, n)
    return(posterior)
  }
    
  
  
  for(i in 1:(k-1)){
    
    n_clust <- rho[i]
    posterior <- posterior + log(theta + i * sigma) + log_pochhammer(1 - sigma, n_clust - 1)
    
    
  }
  
  
  posterior <- posterior + log_pochhammer(1 - sigma, rho[k])
  
  return(posterior)
  
  
  
  
}