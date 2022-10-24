#'@param q: probability to perform a split
#'@param y: data, time X features
#'@param trunc: truncation value transition densities
#'@param step_omega: number of MH steps or number of samples in MC integration
#'


library('extraDistr')
library('MASS')



MCMC <- function(n_iter, burnin, y, q,
                 #method = 'MC_integration',
                 trunc,
                 iter_omega, burnin_omega,
                 alpha_omega, beta_omega, 
                 alpha_sigma, beta_sigma, 
                 alpha_theta, beta_theta){
  

  #### INITIALIZATION ####
  
  prop_sigma <- array(0, n_iter)
  acc_sigma <- 0
  Acc_sigma <- array(0, n_iter)
  
  prop_theta <- array(0, n_iter)
  acc_theta <- 0
  Acc_theta <- array(0, n_iter)
  
  
  ### TODO: save accepted partition
  prop_partition <- list()
  acc_partition_split <- 0
  acc_partition_merge <- 0
  acc_partition_shuffle <- 0
  tot_partition_split <- 0
  tot_partition_merge <- 0
  tot_partition_shuffle <- 0
  Acc_partition <- array(0, n_iter)
  
  
  n <- dim(y)[1]
  d <- dim(y)[2]
  
  # initial partition
  rho <- c(as.integer(n/2))
  rho <- c(rho, n-rho)
  
  # number of clusters
  k <- length(rho)
  
  # initial omega
  omega <- array(0, c(k, d)) # cluster x component
  ### TODO: initialize omega
  
  # initial cluster likelihood
  likelihood <- array(0, k)
  ### TODO: initialize cluster likelihood
  
  
  

  for(iter in (-burnin+1):niter){
    
    # TODO: what and how to save
    
    #### MERGE & SPLIT ####
    
     

    if(runif(1) <= (q * ifelse(k < n, 1, 0) * ifelse(k > 1, 1, 0) + ifelse(k == 1, 1, 0)) ){
      # split
      
      # extreme values are not surely generated (not a.s., but surely by implementation)
      # check help of runif() for details
      
      
      # propose split
      output_split <- split(rho)
      
      rho_proposed <- output_split[[1]]
      j <- output_split[[2]]
      
      
      # compute eppf proposed
      eppf_proposed <- log_EPPF(rho_proposed)
      
      # update omega for new clusters
      omega_proposed <- MH_omega_split(burnin_omega, iter_omega, 
                                  y, j, rho, rho_proposed,
                                  omega, sigma_0, 
                                  alpha_omega, beta_omega, trunc)
      
      # compute likelihood proposed
      likelihood_proposed <- full_log_integrated_likelihood_after_split(likelihood, y, 
                                                                        rho_proposed, j, trunc, 
                                                                        omega_proposed)
      
      # compute log MH-alpha
      
      log_ratio <- MC_log_alpha_split(q, n, j,
                                      likelihood, eppf, 
                                      likelihood_porposed, eppf_proposed,
                                      rho, rho_proposed)
      
      # MH step
      
      tot_split <- tot_split + 1 
      
      if(log(runif(1)) <= min(0, log_ratio)){
        
        rho <- rho_proposed
        omega <- omega_proposed
        likelihood <- likelihood_proposed
        eppf <- eppf_proposed
        acc_partition_split <- acc_partition_split + 1
        
      }
      
      if(iter > 0){
        
        ### TODO: save partition
        
        
      }
      
      
      
      
    } else {
      # merge
      
      # propose merge
      output_list <- merge(rho)
      
      rho_proposed <- output_list[1]
      j <- output_list[2]
      
      # compute eppf proposed
      eppf_proposed <- log_EPPF(rho_proposed)
      
      
      # update omega for new clusters
      omega_proposed <- MH_omega_merge(burnin_omega, iter_omega, 
                                       y, j, d, rho_proposed,
                                       omega, sigma_0, 
                                       alpha_omega, beta_omega, trunc)
      
      # update omega
      
      
      # MH step
      
      
      
      
    }
    
    #### SHUFFLE ####
    
    #### UPDATE PARAMETERS ####
      
   
    
    
    
    
    
  }
  
  
  
}