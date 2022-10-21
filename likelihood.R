#### SUPPORT FUNCTIONS ####
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


#### 
integrated_likelihood_split <- function(type, old_likelihood, y, rho, j, N, omega){
  
  y_partition <- split_data_partition(y, rho)
  
  k <- length(rho)
  
  
  
  
  
  
}