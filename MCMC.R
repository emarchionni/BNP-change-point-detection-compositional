#'@param q: probability to perform a split
#'@param y: data, time X features
#'@param trunc: truncation value transition densities
#'@param step_omega: number of MH steps or number of samples in MC integration
#'


library('extraDistr')
library('MASS')
library('copula')
library('pbmcapply')


source('MCMC_split_merge_shuffle.R')
source('MCMC_alpha.R')
source('MCMC_omega.R')
source('likelihood.R')
source('transition_densities.R')
source('eppf.R')
source('support_functions.R')
source('update_parameters/MCMC_sigma.R')
source('update_parameters/sample_theta.R')



MCMC <- function(niter, burnin, y, q,
                 #method = 'MC_integration',
                 trunc,
                 iter_omega, burnin_omega, sigma_proposal_omega,
                 alpha_omega, beta_omega, 
                 alpha_sigma, beta_sigma,
                 alpha_propose_sigma, beta_propose_sigma,
                 alpha_theta, beta_theta){
  

  #### INITIALIZATION ####
  
  
  
  # saving containers
  Acc_sigma <- array(0, n_iter)
  acc_sigma <- 0
  
  Theta <- array(0, n_iter)
  
  Acc_partition_iter <- list()
  acc_partition_split <- 0
  acc_partition_merge <- 0
  acc_partition_shuffle <- 0
  tot_partition_split <- 0
  tot_partition_merge <- 0
  tot_partition_shuffle <- 0
  
  
  n <- dim(y)[1]
  d <- dim(y)[2]
  
  # initial partition
  #rho <- c(as.integer(n/2))
  #rho <- c(rho, n-rho)
  rho <- rep(1,n)
  
  # number of clusters
  k <- length(rho)
  
  
  # initial omega
  omega <- array(1, c(k, d)) # cluster x component
  # rgamma(d * k, alpha_omega, beta_omega)
  
  # initial cluster likelihood
  likelihood <- full_log_integrated_likelihood(y, rho, omega, n_clust, trunc)
  
  # initial sigma
  sigma <- rbeta(1, alpha_sigma, beta_sigma)
  
  # initial theta
  theta <- rsgamma(-sigma, alpha_theta, beta_theta)
  
  # eppf
  eppf <- log_EPPF(rho, theta, sigma)
  
  # progress bar
  # pb <- progressBar(min = -(burnin-1), max = niter, style = "ETA")
  
  # MCMC
  for(iter in (-burnin+1):niter){
    
    # TODO: what and how to save
    
    print(iter)
    
    #### MERGE & SPLIT ####
    
    k <- length(rho)

    if(runif(1) <= (q * ifelse(k < n, 1, 0) * ifelse(k > 1, 1, 0) + ifelse(k == 1, 1, 0)) ){
      # split
      
      # extreme values are not surely generated (not a.s., but surely by implementation)
      # check help of runif() for details
      
      
      # propose split
      output_split <- split(rho)
      
      rho_proposed <- output_split[[1]]
      j <- output_split[[2]]
      
      
      # compute eppf proposed
      eppf_proposed <- log_EPPF(rho_proposed, theta, sigma)
      
      # update omega for new clusters
      omega_proposed <- MH_omega_split(burnin_omega, iter_omega, 
                                  y, j, d, rho, rho_proposed,
                                  omega, sigma_proposal_omega, 
                                  alpha_omega, beta_omega, trunc)
      
      # compute likelihood proposed
      likelihood_proposed <- full_log_integrated_likelihood_after_split(likelihood, y, 
                                                                        rho_proposed, j, trunc, 
                                                                        omega_proposed)
      
      # compute log MH-alpha
      log_ratio <- MC_log_alpha_split(q, j,
                                      likelihood, eppf, 
                                      likelihood_proposed, eppf_proposed,
                                      rho, rho_proposed)
      
      # MH step
      
      tot_partition_split <- tot_partition_split + 1 
      
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
      
      rho_proposed <- output_list[[1]]
      j <- output_list[[2]]
      
      # compute eppf proposed
      eppf_proposed <- log_EPPF(rho_proposed, theta, sigma)
      
      
      # update omega for new clusters
      omega_proposed <- MH_omega_merge(burnin_omega, iter_omega, 
                                       y, j, rho_proposed,
                                       omega, sigma_proposal_omega, 
                                       alpha_omega, beta_omega, trunc)
      
      
      # compute likelihood proposed
      likelihood_proposed <- full_log_integrated_likelihood_after_merge(likelihood, y, 
                                                                        rho_proposed, j, 
                                                                        trunc, omega_proposed)
      
      
      
      # compute log MH-alpha
      
      
      log_ratio <- MC_log_alpha_merge(q, j,
                                      likelihood, eppf,
                                      likelihood_proposed, eppf_proposed,
                                      rho, rho_proposed)
      
      
      # MH step
      
      tot_partition_merge <- tot_partition_merge + 1
      
      if(log(runif(1)) <= min(0, log_ratio)){
        
        rho <- rho_proposed
        omega <- omega_proposed
        likelihood <- likelihood_proposed
        eppf <- eppf_proposed
        acc_partition_merge <- acc_partition_merge + 1
        
      }
      
      if(iter > 0){
        
        ### TODO: save 
        
      }
      
      
    }
    
    #### SHUFFLE ####
    
    
    if(length(rho) > 1){
      
      
      # propose shuffle
      output_list <- shuffle(rho)
      
      rho_proposed <- output_list[[1]]
      j <- output_list[[2]]
      
      
      
      if(rho_proposed[j] != rho[j]){ 
        # checking whether an actual shuffle occurred
        # to make the code general each function here called checks this condition by its own, 
        # but to fasten it, we check once and for all here
        
        
        # compute eppf proposed
        eppf_proposed <- log_EPPF(rho_proposed, theta, sigma)
        
        
        # update omega for new clusters
        omega_proposed <- MH_omega_shuffle(burnin_omega, iter_omega, 
                                           y, j, d, 
                                           rho, rho_proposed,
                                           omega, sigma_proposal_omega, 
                                           alpha_omega, beta_omega, trunc)
        
        
        # compute likelihood proposed
        likelihood_proposed <- full_log_integrated_likelihood_after_shuffle(likelihood, y, 
                                                                            rho, rho_proposed, j, 
                                                                            trunc, omega_proposed)
        
        
        
        # compute log MH-alpha
        
        
        log_ratio <- MC_log_alpha_shuffle(q, j,
                                          likelihood, eppf,
                                          likelihood_proposed, eppf_proposed,
                                          rho, rho_proposed)
        
        
        # MH step
        
        tot_partition_shuffle <- tot_partition_shuffle + 1
        
        if(log(runif(1)) <= min(0, log_ratio)){
          
          rho <- rho_proposed
          omega <- omega_proposed
          likelihood <- likelihood_proposed
          eppf <- eppf_proposed
          acc_partition_merge <- acc_partition_merge + 1
          
        }
        
        
        if(iter > 0){
          
          #TODO: save
          
        }
        
        
        
      }
    
      
      
      
      
      
    }
    
    
    # save partition
    
    if(iter > 0)
      Acc_partition_iter[[iter]] <- rho
    
    
    #### UPDATE PARAMETERS ####
      
    
    ### sigma
    
    # propose sigma
    sigma_proposed <- rbeta(1, alpha_propose_sigma, beta_propose_sigma)
    
    # compute log MH-alpha
    log_ratio <- MC_log_alpha_sigma(sigma, sigma_proposed,
                                    alpha_sigma, beta_sigma,
                                    theta, alpha_theta, beta_theta,
                                    rho)
    
    if(log(runif(1)) <= min(0, log_ratio)){
      
      sigma <- sigma_proposed
      acc_sigma <- acc_sigma + 1
      
    }
    
    if(iter > 0)
      Acc_sigma[iter] <- sigma
    
      
      
    
    ### theta
    
    # sample theta
    theta <- sample_theta(n, k,
                          theta, sigma, 
                          alpha_theta, beta_theta)
    
    if(iter > 0)
      Theta[iter] <- theta
    
   
    
    # progress bar
    # setTxtProgressBar(pb, iter)
    
  }
  
  
  ### TODO: return output
  
  return(Acc_partition_iter)
}