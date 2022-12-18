### CALIBRATION MC integration ###

# 3-dim simplex


remove(list = ls())

library('extraDistr')
library('MASS')
library('copula')

setwd('C:/Users/edoar/Desktop/Tesi/Code/BNP-change-point-detection-compositional')
source('C:/Users/edoar/Desktop/Tesi/Code/BNP-change-point-detection-compositional/data/death_process.R')




# grid of hyperparams to test
hyperparams <- expand.grid(c(1, 2, 3), c(20, 50, 100, 250, 500)) 
colnames(hyperparams) <- c('trunc', 'iter_omega')
write.table(hyperparams, file = 'simulations/MH_calibration/hyperparameters.txt')
n_couples_param <- dim(hyperparams)[1]

# hyperparams to generate data
n_rep <- 5
n_change_points <- 3
n_data <- 50 # number of observation for each regime

omega_all <- rbind(c(1.3, 0.5, 1),
                   c(0.4, 0.8, 1),
                   c(1.2, 1.5, 0.9))
colnames(omega_all) <- c('1', '2', '3')
write.table(omega_all, file = 'simulations/MH_calibration/omega.txt')
omega_all/apply(omega_all, 1, sum)
trunc <- 50
t <- 1
d <- dim(omega_all)[2]


# non-varying hyperparams of the algorithm
n_iter <- 3000
burnin <- 1000
q <- .5
sigma_proposal_omega <- 0.01
alpha_omega <- 2 
beta_omega <- 2
alpha_sigma <- 1
beta_sigma <- 1
alpha_propose_sigma <- 1 
beta_propose_sigma <- 1
alpha_theta <- 1
beta_theta <- 1


# output list
output_list <- list()


source('MCMC.R')



# loop over replicas
for(rep in 1:n_rep){
  
  # data structure
  y <- array(0, c(0, d))
  
  # GENERATE DATA FOR THE CURRENT REP
  
  for(point in 1:n_change_points){
    
    y_regime <- array(0, c(n_data, d))
    
    omega <- omega_all[point,]
    omega_norm <- sum(omega)
    
    y_old <- rdirichlet(1, omega)
    
    y_regime[1, ] <- y_old
    
    
    for(datum in 2:n_data){
      
      m <- simulate_death_process(trunc, t, omega_norm)
      print(m)
      
      if(m == 0){
        
        l_vect <- rep(0, d) 
        
      } else {
        
        l_vect <- rmnom(1, m, y_old)
        
      }
      
      omega_sim <- omega + l_vect
      
      y_new <- rdirichlet(1, omega_sim)
      
      y_regime[datum, ] <- y_new
      
      y_old <- y_new
      
      
    }
    
    y <- rbind(y, y_regime)
    
  }
  
 
  
  for(idx in 1:n_couples_param){
    
    print(paste('Replica', rep))
    print('Hyperparameters')
    print(paste('trunc:', hyperparams[idx, 1], 'nodes:', hyperparams[idx, 2]))
    
    trunc <- hyperparams[idx, 1]
    iter_omega <- hyperparams[idx, 2]
    
    print(trunc)
    print(iter_omega)
    
    
    
    output_list[[as.character(idx)]][[as.character(rep)]] <- MCMC(n_iter, burnin, y, q,
                                                                  trunc,
                                                                  iter_omega, sigma_proposal_omega,
                                                                  alpha_omega, beta_omega, 
                                                                  alpha_sigma, beta_sigma,
                                                                  alpha_propose_sigma, beta_propose_sigma,
                                                                  alpha_theta, beta_theta, 
                                                                  FALSE, FALSE)
    
    
    
  }
  
  
  
}


save(output_list, file = 'simulations/MC_calibration/output.RData')














