#### POCHHAMMER ####

log_pochhammer <- function(x, factor){
  
  value <- 0
  
  if(factor == 0)
    return(value)
  
  
  for(i  in 1:factor)
    value <- value + log(x + i - 1)
  
  return(value)
  
  
}

pochhammer <- function(x, factor){
  
  value <- 1
  
  
  if(factor == 0)
    return(value)
  
  
  for(i  in 1:factor)
    value <- value * (x + i - 1)
  
  return(value)
  
  
}


#### STIRLING NUMBER ####


abs_stirling_number_first_BC <- function(k, r){
  
  if(k == 0 && r == 0)
    return(1)
  
  if(r <= 0 ||  r > k)
    return(0)
  
  return(abs(Stirling1(k, r)))
  
}


#### SHIFTED GAMMA ####

rsgamma <- function(shifting, alpha, beta){
  
  # shifting is negative
  
  return(rgamma(1, alpha, beta) - shifting)
  
}


#### SPLIT DATA INTO CLUSTERS ####

split_data_partition <- function(y, rho){
  
  y_partition <- list()
  
  for(i in 1:length(rho)){
    
    if(i == 1){
      
      y_partition[[i]] <- y[1:rho[i], ]
      
    } else {
      
      first_index <- sum(rho[1:(i-1)]) + 1
      last_index <- sum(rho[1:(i-1)]) + rho[i]
      
      y_partition[[i]] <- y[first_index:last_index,]
      
    }
    
  }
  
  return(y_partition)
  
}