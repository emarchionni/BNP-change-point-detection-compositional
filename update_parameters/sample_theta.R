sample_theta <- function(n, k,
                         theta, sigma, 
                         alpha_theta, beta_theta){
  
  
  w <- array(0, k + 2)
  
  z <- rbeta(n = 1, theta + 2, n)
  y <- rexp(n = 1, theta + 1)
  
  for(j in 0:(k + 1)){
    
    value <- log(alpha_theta + j) - j * (log(sigma) + log(beta_theta + y - log(z)))
    
    value <- value + log((n - sigma) * (n + 1 - sigma) * abs_stirling_number_first_BC(k - 1, j)
                         + (2 * n + 1 - 2 * sigma) * sigma * abs_stirling_number_first_BC(k - 1, j - 1)
                         + sigma^2 * abs_stirling_number_first_BC(k - 1, j - 2)
                         )
    
    w[j + 1] <- value
    
  }
  
  
  index_chosen <- sample(0:(k+1), size = 1, prob = exp(w))
  
  return(rsgamma(-sigma, alpha_theta + j, beta_theta + t - log(z)))
  
  
  
}