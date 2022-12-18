remove(list = ls())

library('extraDistr')
library('MASS')
library('copula')

setwd('C:/Users/edoar/Desktop/Tesi/Code/BNP-change-point-detection-compositional')
source('C:/Users/edoar/Desktop/Tesi/Code/BNP-change-point-detection-compositional/data/death_process.R')



source('MCMC.R')

omega_all <- rbind(c(1.3, 0.5, 1),
                   c(0.4, 0.8, 1),
                   c(1.2, 1.5, 0.9))
n_change_points <- 3
trunc <- 50
t <- 1
d <- dim(omega_all)[2]
n_data <-50

# data structure
y <- array(0, c(0, d))

# GENERATE DATA

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


n_iter <- 4000
burnin <- 1000
q <- .5
trunc <- 1
iter_omega <- 150

sigma_proposal_omega <- 0.01
alpha_omega <- 2 
beta_omega <- 2
alpha_sigma <- 1
beta_sigma <- 1
alpha_propose_sigma <- 1 
beta_propose_sigma <- 1
alpha_theta <- 1
beta_theta <- 1

output_list <- MCMC(n_iter, burnin, y, q,
                    trunc,
                    iter_omega, sigma_proposal_omega,
                    alpha_omega, beta_omega, 
                    alpha_sigma, beta_sigma,
                    alpha_propose_sigma, beta_propose_sigma,
                    alpha_theta, beta_theta, 
                    FALSE, FALSE)


