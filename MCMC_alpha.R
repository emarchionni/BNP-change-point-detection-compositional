#### ALPHA SPLIT ####

MC_log_alpha_split <- function(q, j,
                               likelihood_old, eppf_old, 
                               likelihood_porposed, eppf_proposed,
                               rho, rho_proposed){
  
  k_proposed <- length(rho_proposed)
  k_old <- length(rho)
  
  n <- sum(rho)
  
  
  if(k_old == 1){
    
    log_ratio <- log(1 - q) + log(n - 1) + sum(likelihood_porposed) + eppf_proposed
    
    log_ratio <- log_ratio - sum(likelihood_old) - eppf_old
    
    return(log_ratio)
    
  }
  
  
  ng <- length(which(rho > 1)) 
  nl <- rho[j]
  
  log_ratio <- log(ng) + log(nl - 1) + log(1 - q) + sum(likelihood_porposed) + eppf_proposed
  
  log_ratio <- log_ratio - log(k_old) - log(q) - sum(likelihood_old) - eppf_old
  
  return(log_ratio)

  
}



#### ALPHA MERGE ####

MC_log_alpha_merge <- function(q, j,
                               likelihood_old, eppf_old,
                               likelihood_proposed, eppf_proposed,
                               rho, rho_proposed){
  
  
  k_proposed <- length(rho_proposed)
  k_old <- length(rho)
  
  n <- sum(rho)
  
  if(k_old == n){
    
    log_ratio <- log(q) + log(n - 1) + sum(likelihood_proposed) + eppf_proposed
    
    log_ratio <- log_ratio - sum(likelihood_old) - eppf_old
    
  }
  
  ng <- length(which(rho_proposed > 1))
  n_new <- rho_proposed[j]
  
  log_ratio <- sum(likelihood_proposed) + eppf_proposed + log(k_old - 1) + log(q)
  
  log_ratio <- log_ratio - sum(likelihood_old) - eppf_old - log(ng) - log(n_new - 1) - log(1 - q)
  
  return(log_ratio)
  
}



#### ALPHA SHUFFLE ####
MC_log_alpha_shuffle <- function(q, j,
                               likelihood_old, eppf_old, 
                               likelihood_porposed, eppf_proposed,
                               rho, rho_proposed){
  
  
  k_proposed <- length(rho_proposed)
  k_old <- length(rho)
  
  n_shuffle <- rho_proposed[j] + rho_proposed[j + 1]
  
  if(n_shuffle == 2)  # i.e. no actual shuffle occurred
    return(0)
  
  
  log_ratio <- sum(likelihood_porposed) + eppf_proposed - sum(likelihood_old) - eppf_old
  
  return(log_ratio)
  
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  