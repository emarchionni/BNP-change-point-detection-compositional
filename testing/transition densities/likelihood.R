#### LOG INTEGRATED LIKELIHOOD CLUSTER ####

log_integrated_likelihood_cluster_cpp <- function(y, omega, n_clust, trunc){
  
  #'@param y: observation in the cluster / dim: time x components
  #'@param omega: omega of cluster / dim: d
  
  if(n_clust == 1)
    return(ddirichlet(y, omega, log = T))
  
  #value <- ddirichlet(y[1,], omega, log = T)
  #
  #for(i in 2:n_clust){
  #  value <- value + log_transition_densities(y[i-1, ], y[i, ], omega, trunc)
  #  if(is.nan(value)) browser()
  #}
  #return(value)
  
  return(log_integrated_likelihood_cluster_multiple_cpp(y, omega, trunc))

  
}


#### LOG INTEGRATED LIKELIHOOD CLUSTER ####

log_integrated_likelihood_cluster <- function(y, omega, n_clust, trunc){
  
  #'@param y: observation in the cluster / dim: time x components
  #'@param omega: omega of cluster / dim: d
  
  if(n_clust == 1)
    return(ddirichlet(y, omega, log = T))
  
  value <- ddirichlet(y[1,], omega, log = T)
  
  for(i in 2:n_clust){
    value <- value + log_transition_densities(y[i-1, ], y[i, ], omega, trunc)
  }
  
  
  
  
  return(value)
  
}

