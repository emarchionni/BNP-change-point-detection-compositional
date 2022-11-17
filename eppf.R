#### SUPPORT FUNCTIONS ####

log_EPPF <- function(rho, theta, sigma){
  
  n <- sum(rho)
  k <- length(rho)

  value <- lfactorial(n) - lfactorial(k) - log_pochhammer(theta + 1, n - 1)
  
  if(k == 1){
    value <- value + log(theta + sigma) + log_pochhammer(1 - sigma, n - 1) + lfactorial(n)
    return(value)
  }
    
  
  for (j in 1:(k-1)) {
  
    nj <- rho[j]
    
    value <- value + log(theta + j * sigma) + log_pochhammer(1 - sigma, nj - 1) - lfactorial(nj)
    
  }
  
  nj <- rho[k]
  
  value <- value + log_pochhammer(1 - sigma, nj - 1) - lfactorial(nj)
  
  return(value)
  
}