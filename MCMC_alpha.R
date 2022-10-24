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
    
    log_ratio <- log_ratio - sum(likelihood_old) + eppf_old
    
    return(ratio)
    
  }
  
  
  ng <- length(which(rho > 1)) 
  nl <- rho[j]
  
  log_ratio <- log(ng) + log(nl - 1) + log(1 - q) + sum(likelihood_porposed) + eppf_proposed
  
  log_ratio <- log_ratio - log(k_old) - log(q) - sum(likelihood_old) - eppf_old
  
  return(log_ratio)

  
}


