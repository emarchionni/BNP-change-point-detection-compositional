#### LOG INTEGRATED LIKELIHOOD CLUSTER ####
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

log_integrated_likelihood_cluster <- function(y, omega, n_clust, trunc){
  
  #'@param y: observation in the cluster / dim: time x components
  #'@param omega: omega of cluster / dim: d
  
  if(n_clust == 1)
    return(ddirichlet(y,omega, log = T))
  
  value <- ddirichlet(y[1,], omega, log = T)
  
  for(i in 2:n_clust)
    value <- value + log(transition_densities(y[i-1, ], y[i, ], omega, trunc))
  
  return(value)
  
}




#### LOG INTEGRATED LIKELIHOOD SPLIT ####
full_log_integrated_likelihood_after_split <- function(old_likelihood, y, 
                                                       rho_proposed, j, 
                                                       trunc, omega){
  
  #'@param omega: updated omega / dim: new_clusters x components
  #'@param old_likelihood: nonupdated likelihood / dim: old_likelihood
  
  y_partition <- split_data_partition(y, rho_proposed)
  

  
  
  new_value <- array(0, 2)
  
  
  for(clust in 0:1){
    
    n_clust <- rho_proposed[j + clust]
    y_clust <- y_partition[[j + clust]]
    omega_clust <- omega[j + clust, ]
    
    
    new_value[clust + 1] <- log_integrated_likelihood_cluster(y_clust, omega_clust, n_clust, trunc)
    
  }
  
  k <- length(rho_proposed) - 1 # old number of clusters
  
  if(k == 1){
    
    return(new_value)
    
  }
  
  if(j == 1){
    
    new_likelihood <- c(new_value, old_likelihood[2:k])
    
  } else if(j == k){
    
    new_likelihood <- c(old_likelihood[1:(k-1)], new_value)
    
  } else {
    
    new_likelihood <- c(old_likelihood[1:(j-1)], new_value, old_likelihood[(j+1):k])
    
  }
  
  if(length(new_likelihood) != length(old_likelihood) + 1)
    browser()
  
  
  return(new_likelihood)
  
  
}


#### LOG INTEGRATED LIKELIHOOD MERGE ####
full_log_integrated_likelihood_after_merge <- function(old_likelihood, y, 
                                                       rho_proposed, j, 
                                                       trunc, omega){
  
  #'@param omega: updated omega / dim: new_clusters x components
  #'@param old_likelihood: nonupdated likelihood / dim: old_likelihood

  y_partition <- split_data_partition(y, rho_proposed)
  
  n_clust <- rho_proposed[j]
  y_clust <- y_partition[[j]]
  omega_clust <- omega[j, ]
  
  new_value <- log_integrated_likelihood_cluster(y_clust, omega_clust, n_clust, trunc)
  
  
  k <- length(rho_proposed) + 1 # old number of clusters
  
  if(k == 2){
    
    return(new_value)
    
  }
  
  if(j == 1){
    
    new_likelihood <- c(new_value, old_likelihood[3:k])
    
  } else if(j == (k-1)){
    
    new_likelihood <- c(old_likelihood[1:(j-1)], new_value)
    
  } else {
    
    new_likelihood <- c(old_likelihood[1:(j-1)], new_value, old_likelihood[(j+2):k])
    
  }
  
  
}


#### LOG INTEGRATED LIKELIHOOD SHUFFLE ####
full_log_integrated_likelihood_after_shuffle <- function(old_likelihood, y, 
                                                         rho_proposed, j, 
                                                         trunc, omega){
  
  
  #'@param omega: updated omega / dim: clusters x components
  #'@param old_likelihood: nonupdated likelihood / dim: likelihood
  
  
  n_shuffle <- rho_proposed[j] + rho_proposed[j + 1]
  
  if(n_shuffle == 2) # i.e. no actual shuffle occurred
    return(old_likelihood)
  
  
  
  y_partition <- split_data_partition(y, rho_proposed)
  
  
  for(clust in 0:1){
    
    n_clust <- rho_proposed[j + clust]
    y_clust <- y_partition[[j + clust]]
    omega_clust <- omega[j + clust, ]
    
    
    old_likelihood[j + clust] <- log_integrated_likelihood_cluster(y_clust, omega_clust, n_clust, trunc)
    
  }
  
  return(old_likelihood)
  
}
