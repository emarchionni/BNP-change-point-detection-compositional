


pochhammer <- function(x, factor){
  
  value <- 1
  
  
  if(factor == 0)
    return(value)
  
  
  for(i  in 1:factor)
    value <- value * (x + i - 1)
  
  return(value)
  
  
}


simulate_death_process <- function(trunc, t, omega_norm){
  
  q <- array(0, trunc + 1)
  
  # for q0
  value <- 0
  for(n in 1:trunc){
    temp <- exp( - 0.5 * t * n * (n - 1 + omega_norm)) * (-1)^n * (omega_norm + 2 * n - 1) * pochhammer(omega_norm, n-1)
    temp <- temp / factorial(n) 
    value <- value + temp
  }
    
  q[1] <- 1 - value
  
    
  for(l in 1:trunc){
    
    value <- 0
    
    for(n in l:trunc){
      
      temp <- exp( - 0.5 * t * n * (n - 1 + omega_norm)) * (-1)^(n - l) * choose(n, l)
      temp <- temp * pochhammer(omega_norm + l, n - 1) * (omega_norm + 2 * n - 1) / factorial(n)
      value <- value + temp
      
    }
    
    q[l + 1] <- value
    
  }
  
  return(sample(0:trunc, size = 1, prob = q))
  
}