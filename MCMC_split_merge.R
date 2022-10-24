#### SPLIT ####


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
    
    rho_proposed <- c(rho[1:(k-1)], l, rho[k] - l)
    
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




split <- function(rho){
  
  k <- length(rho)
  
  
  if(length(rho) == 1){ ### all observation are in a single cluster
    
    return(split_singlecluster(rho))
    
  } else { ### observations are split at least in two clusters
    
    return(split_multicluster(rho))
  }
  
} 
  

  
  
#### MERGE ####

merge_twocluster <- function(rho){
  
  j <- 1
  
  return(list(sum(rho), j))
  
}


merge_multicluster <- function(rho){
  
  
  k <- length(rho)
  
  j <- sample(1:(k-1), size = 1)
  
  
  if(j == 1){
    
    rho_proposed <- c(rho[1] + rho[2], rho[3:k])
    
  } else if (j == (k-1)) {
    
    rho_proposed <- c(rho[1:(j-1)], rho[j] + rho[j+1])
    
  } else {
    
    rho_proposed <- c(rho[1:(j-1)], rho[j] + rho[j+1], rho[(j+2):k])
    
  }
  
  return(list(rho_proposed, j))
  
}


merge <- function(rho){
  
  k <- length(rho)
  
  if(k == 2){ ### observations are split in TWO clusters
    
    return(merge_twocluster(rho))
    
  } else { ### observations are split at least in more than two clusters
    
    return(merge_multicluster(rho))
    
    }
 
  
}
  

