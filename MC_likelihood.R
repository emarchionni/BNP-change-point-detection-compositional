#### LOG INTEGRATED MARGINAL LIKELIHOOD CLUSTER ####

MC_log_integrated_MARGINAL_likelihood_cluster <- function(y, n_clust, trunc, 
                                                          iter_omega,
                                                          alpha_omega, beta_omega){
  
  #'@param y: observation in the cluster / dim: time x components
  
  if(n_clust == 1)
    return(MC_log_integrated_MARGINAL_likelihood_cluster_single_cpp(y, trunc, 
                                                                    iter_omega,
                                                                    alpha_omega, beta_omega))
  
  #value <- ddirichlet(y[1,], omega, log = T)
  #
  #for(i in 2:n_clust){
  #  value <- value + log_transition_densities(y[i-1, ], y[i, ], omega, trunc)
  #  if(is.nan(value)) browser()
  #}
  #return(value)
  
  return(MC_log_integrated_MARGINAL_likelihood_cluster_multiple_cpp(y, trunc, 
                                                                    iter_omega,
                                                                    alpha_omega, beta_omega))
  
  
}


#### INITIALIZATION MARGINAL LIKELIHOOD ####

MC_full_log_integrated_MARGINAL_likelihood <- function(y, rho, trunc, n_clust, 
                                                       iter_omega,
                                                       alpha_omega, beta_omega){
  
  y_partition <- split_data_partition(y, rho)
  
  k <- length(rho)
  
  likelihood <- array(0, k)
  
  for (i in 1:k) {
    
    y_clust <- y_partition[[i]]
    n_clust <- rho[i]

    
    likelihood[i] <- MC_log_integrated_MARGINAL_likelihood_cluster(y_clust, n_clust, trunc, 
                                                                   iter_omega,
                                                                   alpha_omega, beta_omega)
    
  }
  
  return(likelihood)
  
}



#### LOG INTEGRATED MARGINAL LIKELIHOOD SPLIT ####
MC_full_log_integrated_MARGINAL_likelihood_after_split <- function(old_likelihood, y,
                                                                   rho_proposed, j,
                                                                   trunc,
                                                                   iter_omega,
                                                                   alpha_omega, beta_omega){
  
  
  #'@param old_likelihood: nonupdated likelihood / dim: old_likelihood
  
  y_partition <- split_data_partition(y, rho_proposed)
  
  new_value <- array(0, 2)

  for(clust in 0:1){
    
    n_clust <- rho_proposed[j + clust]
    y_clust <- y_partition[[j + clust]]

    
    new_value[clust + 1] <- MC_log_integrated_MARGINAL_likelihood_cluster(y_clust, n_clust, trunc, 
                                                                          iter_omega,
                                                                          alpha_omega, beta_omega);
    
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


#### LOG INTEGRATED MARGINAL LIKELIHOOD MERGE ####
MC_full_log_integrated_MARGINAL_likelihood_after_merge <- function(old_likelihood, y, 
                                                                   rho_proposed, j, 
                                                                   trunc,
                                                                   iter_omega,
                                                                   alpha_omega, beta_omega){
  
  #'@param omega: updated omega / dim: new_clusters x components
  #'@param old_likelihood: nonupdated likelihood / dim: old_likelihood
  
  k_proposed <- length(rho_proposed)
  
  y_partition <- split_data_partition(y, rho_proposed)
  
  n_clust <- rho_proposed[j]
  y_clust <- y_partition[[j]]
  

  
  
  new_value <- MC_log_integrated_MARGINAL_likelihood_cluster(y_clust, n_clust, trunc, 
                                                             iter_omega,
                                                             alpha_omega, beta_omega);
  
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


#### LOG INTEGRATED MARGINAL LIKELIHOOD SHUFFLE ####
MC_full_log_integrated_MARGINAL_likelihood_after_shuffle <- function(old_likelihood, y,
                                                                     rho, rho_proposed, j, 
                                                                     trunc, 
                                                                     iter_omega,
                                                                     alpha_omega, beta_omega){
  
  
  #'@param omega: updated omega / dim: clusters x components
  #'@param old_likelihood: nonupdated likelihood / dim: likelihood
  
  
  
  
  if(rho[j] == rho_proposed[j]) # i.e. no actual shuffle occurred
    return(old_likelihood)
  
  
  
  y_partition <- split_data_partition(y, rho_proposed)
  
  
  for(clust in 0:1){
    
    n_clust <- rho_proposed[j + clust]
    y_clust <- y_partition[[j + clust]]

    
    
    old_likelihood[j + clust] <- MC_log_integrated_MARGINAL_likelihood_cluster(y_clust, n_clust, trunc, 
                                                                               iter_omega,
                                                                               alpha_omega, beta_omega)
    
  }
  
  return(old_likelihood)
  
}
