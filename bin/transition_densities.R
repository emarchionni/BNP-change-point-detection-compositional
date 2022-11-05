#### SUPPORT FUNCTIONS ####



# Zeta_m

zeta_m <- function(y_0, y, omega, m){
  
  if(m == 0)
    return(1)
  
  d <- length(y)
  
  composition <- get_all_weak_composition(m, d)
  
  if(choose(m+d-1, m) == 1) # it should not happen, if so next line does not work -> at least d compositions
    browser()
    
  

  num_composition <- dim(composition)[1]
  
  value <- 0
  
  for (i in 1:num_composition) {
    
    curr_composition <- composition[i,]
    
    value <- value + factorial(m) * (sum(y_0))^curr_composition[d] * ddirichlet(y, omega + curr_composition) / prod(factorial(curr_composition))
    
  }
  
  
  return(value)
  
  
}


# Q_n
Q_n_poly <- function(y_0, y, omega, n){
  
  
  if(n == 0) # Q_0 = 1
    return(1)
  
  
  omega_norm <- sum(omega)
  value <- 0
  
  
  for(m in 0:n){
    
    p <- pochhammer(omega_norm + m, n - 1)
    z_m <- zeta_m(y_0, y, omega, m)
    
  
    value <- value + ((-1)^(n - m)) * choose(n, m) * p * z_m
    
    
  }
  
  
  
  return((omega_norm + 2 * n - 1) * value / factorial(n))
    
}


#### TRANSITION DENSITIES FUNCTION ####

log_transition_densities <- function(y_0, y, omega, trunc){
  
  value_1 <- ddirichlet(y, omega)
  Qn <- array(0, trunc + 1)
  lambda_n <- array(0, trunc + 1)
  
  omega_norm <- sum(omega)
  
  for(n in 0:trunc){
    
    Qn[n + 1] <- Q_n_poly(y_0, y, omega, n)
    lambda_n[n + 1] <- (.5) * n * (n - 1 + omega_norm)
    
  }
  
  
  ## TODO: how to fix this? 
  if((sum(exp(- lambda_n) * Qn) + value_1) < 0){
    warning('transition 0 probability')
    browser()
  }
  
  
  return(log(sum(exp(- lambda_n) * Qn) + value_1))
  
  
}
