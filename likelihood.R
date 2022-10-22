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


#### SPLIT ####
log_integrated_likelihood_after_split <- function(type, old_likelihood, y, rho, j, trunc, omega){
  
  y_partition <- split_data_partition(y, rho)
  k <- length(rho)
  
  
  new_value <- array(0, 2)
  
  
  for(clust in 0:1){
    
    if(j == k)
      j <- k-1
    
    y_cluster <- y_partition[[j + clust]]
    
    
    if(rho[j + clust] == 1){
      
      new_value[clust + 1] <- ddirichlet(y_cluster, omega, log = T)
      
    } else {
      
      new_value[clust + 1] <- ddirichlet(y_cluster[1, ], omega, log = T)
      
      n_clust <- dim(y_cluster)[1]
      
      for(i in 2:n_clust)
        new_value[clust + 1] <- new_value[clust + 1] + log(transition_densities(y_cluster[i-1,], y_cluster[i,], omega, trunc))
      
    }
    
  }
  
  
  if(j == 1){
    
    new_likelihood <- c(new_value, old_likelihood[2:k])
    
  } else if(j == k){
    
    new_lkelihood <- c(old_likelihood[1:(k-1)], new_value)
    
  } else {
    
    new_likelihood <- c(old_likelihood[1:(j-1)], new_value, old_likelihood[(j+1):k])
    
  }
  
  if(length(new_likelihood) != length(old_likelihood) + 1)
    browser()
  
  return(new_likelihood)
  
  
}