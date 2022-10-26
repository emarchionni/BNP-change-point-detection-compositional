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

#### ALPHA THETA ####