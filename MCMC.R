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
  
  prop_partition <- list()
  acc_partition <- 0
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
      
      # update omega for new clusters
      
      
      # compute MH alpha
      
      
      
      # MH step
       
      
      
      
      
    } else {
      # merge
      
      # propose merge
      # update omega
      # MH step
      
    }
    
    #### SHUFFLE ####
    
    #### UPDATE PARAMETERS ####
      
   
    
    
    
    
    
  }
  
  
  
}