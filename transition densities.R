#### SUPPORT FUNCTIONS ####


# weak compositions

get_first_weak_composition <- function(m, d)
{
  composition <- array(0,c(1,d))
  
  for (i in 1:(d-1)) {
    
    composition[i] = 0
    
  }
  composition[d] = m
  #browser()
  return(composition)
}


get_next_weak_composition <- function(m, d, composition)
{
  if (composition[1] == m){
    
    return(F)
    
  }
  
  # there's an i with composition[i] > 0, and it is not 0. find the last one
  last = d
  
  while (composition[last] == 0) {
    
    #browser()
    last <- last-1
    
  }
  
  
  z = composition[last]
  composition[last - 1] = composition[last - 1] + 1
  composition[last] = 0
  composition[d] = z - 1
  return(composition)
  
}



get_all_weak_composition <- function(m, d) {
  
  n_vec <- array(0, c(0, d))
  
  
  composition <- get_first_weak_composition(m, d)
  
  while(is.array(composition)){
    
    if(is.array(composition)){
      
      n_vec <- rbind(n_vec, composition)
      
    }
    
    composition <- get_next_weak_composition(m, d, composition)
    
  }
  
  
  
  if(dim(n_vec)[1]!=choose(m+d-1, m))
    warning('Dimension of array and number of weak compositions do not agree')
  
  return(n_vec)
  
}

# Pochhammer symbol

pochhammer <- function(x, factor){
  
  value <- 1
  
  if(factor == 0)
    return(value)
  
  
  for(i  in 1:factor)
    value <- value * (x + i - 1)
  
  return(value)
  
  
}

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
transition_densities <- function(y_0, y, omega, trunc){
  
  value_1 <- ddirichlet(y, omega)
  log_Qn <- array(0, trunc + 1)
  lambda_n <- array(0, trunc + 1)
  
  omega_norm <- sum(omega)
  
  for(n in 1:trunc){
    log_Qn[n + 1] <- log(Q_n_poly(y_0, y, omega, n))
    lambda_n[n + 1] <- (.5) * n * (n - 1 + omega_norm)
  }
    
  
    
  # browser()
  return(sum(exp(-lambda_n + log_Qn)) + value_1)
    
  
}