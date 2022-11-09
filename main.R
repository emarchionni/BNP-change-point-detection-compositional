#### SIMULATED DATA ####

remove(list = ls())


library('extraDistr')
library('MASS')
library('copula')

omega1 <- c(2, 5, 10)
omega1_mean <- omega1/sum(omega1)
y1 <- rdirichlet(50, omega1)

omega2 <- c(10, 5, 5)
omega2_mean <- omega2/sum(omega2)
y2 <- rdirichlet(50, omega2)

omega3 <- c(3, 1, .2)
omega3_mean <- omega3/sum(omega3)
y3 <- rdirichlet(50, omega3)

omega1_mean; omega2_mean; omega3_mean;

y <- rbind(y1, y2, y3)

remove(list = c('omega1', 'y1', 'omega2', 'y2', 'omega3', 'y3'))

setwd('C:/Users/edoar/Desktop/Tesi/Code/BNP-change-point-detection-compositional')

n_iter <- 1000
burnin <- 0
q <- .8
trunc <- 5
iter_omega <- 100
#burnin_omega <- 50
sigma_proposal_omega <- 0.01
alpha_omega <- 2 
beta_omega <- 2
alpha_sigma <- 1
beta_sigma <- 1
alpha_propose_sigma <- 1 
beta_propose_sigma <- 1
alpha_theta <- 1
beta_theta <- 1

source('MCMC.R')



output <- MCMC(n_iter, burnin, y, q,
                 #method = 'MC_integration',
                 trunc,
                 iter_omega, sigma_proposal_omega,
                 alpha_omega, beta_omega, 
                 alpha_sigma, beta_sigma,
                 alpha_propose_sigma, beta_propose_sigma,
                 alpha_theta, beta_theta)


