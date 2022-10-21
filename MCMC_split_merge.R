#### SPLIT ####

split <- function(rho){
  
  k <- length(rho)
  
  
  if(length(rho) == 1){ ### all observation were in a single cluster
    
    return(split_singlecluster(rho))
    
  } else { ### observations were split at least in two clusters
    
    return(split_multicluster(rho))
  }
  
} 
  

split_singlecluster <- function(rho){
  
  # choose where to split the unique cluster (selecting the dim of the new cluster)
  l <- sample(1:(rho[1]-1), size = 1)
  
  rho_proposed <- c(l, rho[1]-l)
  
  return(list(rho_proposed,1))
  
}
  

split_multicluster <- function(rho){
  
  k <- length(rho)
  
  # select a cluster s.t. n_j > 1 
  idx_no_single <- which(rho > 1)
  if(length(idx_no_single) != 1){
    
    j <- sample(idx_no_single, size = 1)
    
  } else {
    
    j <- idx_no_single
    
  }
  
  
  # choose where to split the picked cluster (selecting the dim of the new cluster)
  l <- sample(1:(rho[j]-1), size = 1)
  
  # propose new partition
  if(j == 1){
    
    rho_proposed <- c(l, rho[1] - l, rho[2:k])
    
  } else if(j == k) {
    
    rho_proposed <- c(rho[1:k-1], l, rho[k] - l)
    
  } else {
    
    rho_proposed <- c(rho[1:(j-1)], l, rho[j] - l, rho[(j+1):k])
    
  }
  
  # check code
  if(length(rho_proposed) != k+1){
  
    warning('Problem partition proposed')
    browser()
    
  }
  
  
  return(list(rho_proposed,j))
  
  
}  
  
  
  
  
  

