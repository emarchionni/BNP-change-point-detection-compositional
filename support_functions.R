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


#### GET WEAK COMPOSITIIONS OF INTEGERS ####

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



