pochhammer <- function(x, factor){
  
  value <- 1
  
  
  if(factor == 0)
    return(value)
  
  
  for(i  in 1:factor)
    value <- value * (x + i - 1)
  
  return(value)
  
  
}



b_k_m <- function(k, m, theta, t){
  
  value <- pochhammer(theta + m, k - 1) * (theta + 2 * k - 1) / factorial(m)
  
  return( exp(- 0.5 * k * (k - 1 + theta) * t) * value )
  
}



S_minus_k <- function(k_vect, theta, t){
  
  value <- 0
  m <- length(k_vect) 
  
  for(i in 1:m)
    for(j in 1:(2 * k_vect[i] + 1))
      value <- value + (-1)^j * b_k_m(m + i, m, theta, t)
      
    
  return(value)
  
  
}


S_plus_k <- function(k_vect, theta, t){
  
  value <- 0
  m <- length(k_vect)
  
  for(i in 1:m)
    for(j in 1:(2 * k_vect[i]))
      value <- value + (-1)^j * b_k_m(m + i, m, theta, t)
  
  
  return(value)
  
  
}



simulating_death_process <- function(omega_norm, t){
  
  m <- 0
  k_m <- 0
  k_vect <- c(k_m)
  
  u <- runif(1)
  
  while(TRUE){
    
    i <- 0
    
    while(b_k_m(i + m + 1, m, omega_norm, t) < b_k_m(i + m, m, omega_norm, t))
      i <- i + 1
    
    k_m <- ceiling(i * 0.5)
    

    
    k_vect_computation <- k_vect
    
    while((S_minus_k(k_vect_computation, omega_norm, t) < u) && (u < S_plus_k(k_vect_computation, omega_norm, t)))
      k_vect_computation <- k_vect_computation + rep(1, length(k_vect_computation))
    
    if(S_minus_k(k_vect_computation, omega_norm, t) > u){
      
      return(m)
      
    } else {
      
      if(m == 0)
        k_vect <- c(k_vect, k_m)
      else
        k_vect <- c(k_vect[1:(length(k_vect)-1)], k_m, 0)
  
      m <- m+1
      
      
    }
    
    
    
    
    
    
    
  }
  
  
  
  
}
