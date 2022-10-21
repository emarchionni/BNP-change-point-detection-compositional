#### ALPHA ####

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


MC_alpha <- function(trunc, step_omega, y, omega, j, rho_proposed, rho){
  
  k_proposed <- length(rho_proposed)
  k_previous <- length(rho)
  
  y_proposed <- split_data_partition(y, rho_proposed)
  y_previous <- split_data_partition(y, rho)
    
    

  
}


