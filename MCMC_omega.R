#### MH SINGLE CLUSTER OMEGA ####



MH_omega <- function(burnin_omega, iter_omega, 
                     y, 
                     initial_omega, sigma_0, 
                     alpha_omega, beta_omega, 
                     n_clust, trunc){
  
  #'@param y: observation in the cluster / dim: time x components
  #'@param omega_initial: old_omega of cluster, where to centre gaussian / dim: d
  
  d <- length(initial_omega)
  
  sample_omega <- array(0, c(iter_omega, d))
  
  sigma_matrix <- diag(d) * sigma_0
  
  old_omega <- log(initial_omega)
  
  for(i in (-burnin_omega+1):iter_omega){
    
    proposed_omega <- mvrnorm(n = 1, log(initial_omega), sigma_matrix)
    
    numerator <- log_integrated_likelihood_cluster(y, exp(proposed_omega), n_clust, trunc)
    denominator <- log_integrated_likelihood_cluster(y, exp(old_omega), n_clust, trunc)
    
    for(j in 1:d){
      if(numerator == Inf || denominator == Inf) browser()
      numerator <- numerator + dgamma(exp(proposed_omega[j]), shape = alpha_omega, rate = beta_omega, log = T)
      denominator <- denominator + dgamma(exp(old_omega[j]), shape = alpha_omega, rate = beta_omega, log = T)
      
    }
    
    # det of Jacobian
    numerator <- numerator + sum(proposed_omega)
    denominator <- denominator + sum(old_omega)
    
    ratio <- numerator - denominator
    
    #if(is.nan(ratio)) browser()
    
    if(log(runif(1)) <= min(0, ratio)){
      old_omega <- proposed_omega
      #print('accepted')
    } #else print('refused')
      
    
    if(i > 0)
      sample_omega[i, ] <- old_omega
      
    
  }
  
  
  return(colMeans(exp(sample_omega)))
  
}



#### MH OMEGA SPLIT ####

MH_omega_split <- function(burnin_omega, iter_omega, 
                           y, j, d, rho, rho_proposed,
                           omega, sigma_0, 
                           alpha_omega, beta_omega, trunc){
  
  #'@param omega: nonupdated omega / dim: old_clusters x components
  #'@param j: j and j+1 are the new split clusters
  
  y_partition <- split_data_partition(y, rho_proposed)

  
  new_value <- array(0, c(2,d))
  
  
  k_old <- length(rho)
  
  if(k_old == 1){
    
    initial_omega <- omega
    
  } else {
    
    initial_omega <- omega[j, ]
    
  }
  
  
  
  for(clust in 0:1){
    
    n_clust <- rho_proposed[j + clust]
    y_clust <- y_partition[[j + clust]]
    
    
    new_value[clust + 1, ] <- MH_omega(burnin_omega, iter_omega, y_clust, 
                                       initial_omega, sigma_0, 
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

MH_omega_merge <- function(burnin_omega, iter_omega, 
                           y, j, rho_proposed,
                           omega, sigma_0, 
                           alpha_omega, beta_omega, trunc){
  
  #'@param omega: nonupdated omega / dim: old_clusters x components
  #'@param j: merged cluster
  
  #d <- dim(omega)[2] #there were at least two clusters
  
  
  y_partition <- split_data_partition(y, rho_proposed)
  
  initial_omega <- colMeans(omega[j:(j+1), ]) # initialize omega on the mean of the omega of the two merged clusters
  
  n_clust <- rho_proposed[j]
  y_clust <- y_partition[[j]]
  
  
 new_value <-  MH_omega(burnin_omega, iter_omega, y_clust, 
                        initial_omega, sigma_0, 
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

MH_omega_shuffle <- function(burnin_omega, iter_omega, 
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
 
    
    omega[j + clust, ] <- MH_omega(burnin_omega, iter_omega, y_clust, 
                                   omega[j + clust, ], sigma_0, 
                                   alpha_omega, beta_omega, 
                                   n_clust, trunc)
    
  }
  
  
  
  return(omega)
  
  
  
  
}
