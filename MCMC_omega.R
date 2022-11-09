#### MH SINGLE CLUSTER OMEGA ####



MH_omega <- function(iter_omega, 
                     y, 
                     d, sigma_0, 
                     alpha_omega, beta_omega, 
                     n_clust, trunc){
  
  #'@param y: observation in the cluster / dim: time x components
  
  if(n_clust == 1)
    return(MH_omega_cluster_single(iter_omega, y, d, sigma_0, alpha_omega, beta_omega, trunc))
  
  
  return(MH_omega_cluster_multi(iter_omega, y, d, sigma_0, alpha_omega, beta_omega, trunc))
  
  
}



#### MH OMEGA SPLIT ####

MH_omega_split <- function(iter_omega, 
                           y, j, d, rho_proposed,
                           omega, sigma_0, 
                           alpha_omega, beta_omega, trunc){
  
  #'@param omega: nonupdated omega / dim: old_clusters x components
  #'@param j: j and j+1 are the new split clusters
  
  y_partition <- split_data_partition(y, rho_proposed)

  
  new_value <- array(0, c(2,d))
  
  
  
  
  for(clust in 0:1){
    
    n_clust <- rho_proposed[j + clust]
    y_clust <- y_partition[[j + clust]]
    
    
    new_value[clust + 1, ] <- MH_omega(iter_omega, y_clust, 
                                       d, sigma_0, 
                                       alpha_omega, beta_omega, 
                                       n_clust, trunc)
     
    
  }
  
  k <- length(rho_proposed) - 1 # old number of clusters
  
  if(k == 1){
    
    return(new_value)
    
  }
  
  if(j == 1){
    
    new_omega <- rbind(new_value, omega[2:k,])
    
  } else if(j == k){
    
    new_omega <- rbind(omega[1:(k-1),], new_value)
    
  } else {
    
    new_omega <- rbind(omega[1:(j-1),], new_value, omega[(j+1):k,])
    
  }
  

  
  
  return(new_omega)
  
  
  
}


#### MH OMEGA MERGE ####

MH_omega_merge <- function(iter_omega, 
                           y, j, rho_proposed,
                           omega, sigma_0, 
                           alpha_omega, beta_omega, trunc){
  
  #'@param omega: nonupdated omega / dim: old_clusters x components
  #'@param j: merged cluster
  
  d <- dim(omega)[2] #there were at least two clusters
  
  
  y_partition <- split_data_partition(y, rho_proposed)
  
  
  n_clust <- rho_proposed[j]
  y_clust <- y_partition[[j]]
  
  
 new_value <-  MH_omega(iter_omega, y_clust, 
                        d, sigma_0, 
                        alpha_omega, beta_omega, 
                        n_clust, trunc)
 
 
 k <- length(rho_proposed) + 1 # old number of clusters
 
 if(k == 2){
   
   return(new_value)
   
 }
 
 
 if(j == 1){
   
   new_omega <- rbind(new_value, omega[3:k, ])
   
 } else if (j == (k-1)){
   
   new_omega <- rbind(omega[1:(k-2), ], new_value)
   
 } else {
   
   new_omega <- rbind(omega[1:(j-1), ], new_value, omega[(j+2):k, ])
   
 }
   
  
  return(new_omega)
  
  
}


#### MH OMEGA SHUFFLE ####

MH_omega_shuffle <- function(iter_omega, 
                           y, j, d, 
                           rho, rho_proposed,
                           omega, sigma_0, 
                           alpha_omega, beta_omega, trunc){
  
  #'@param omega: nonupdated omega / dim: old_clusters x components
  #'@param j: j and j+1 are the shuffled clusters
  
  
  if(rho[j] == rho_proposed[j]) # i.e. no actual shuffle occurred
    return(omega)
  
  y_partition <- split_data_partition(y, rho_proposed)
  
  
  new_value <- array(0, c(2,d))
  
  
  for(clust in 0:1){
    
    y_clust <- y_partition[[j + clust]]
    n_clust <- rho_proposed[j + clust]
 
    
    omega[j + clust, ] <- MH_omega(iter_omega, y_clust, 
                                   d, sigma_0, 
                                   alpha_omega, beta_omega, 
                                   n_clust, trunc)
    
  }
  
  
  
  return(omega)
  
  
  
  
}
