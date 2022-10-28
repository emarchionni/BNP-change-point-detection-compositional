#### SUPPORT FUNCTIONS ####

log_EPPF <- function(rho, theta, sigma){
  
  n <- sum(rho)
  k <- length(rho)

  value <- log(factorial(n)) - log(factorial(k)) - log(pochhammer(theta + 1, n - 1))
  
  if(k == 1){
    value <- value + log(theta + sigma) + log(pochhammer(1 - sigma, n - 1)) + log(factorial(n))
    return(value)
  }
    
  
  for (j in 1:(k-1)) {
  
    nj <- rho[j]
    
    value <- value + log(theta + j * sigma) + log(pochhammer(1 - sigma, nj - 1)) - log(factorial(nj))
    
  }
  
  nj <- rho[k]
  
  value <- value + log(pochhammer(1 - sigma, nj - 1)) - log(factorial(nj))
  
  return(value)
  
}