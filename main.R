#### SIMULATED DATA ####

remove(list = ls())


library('extraDistr')
library('MASS')
library('copula')

omega1 <- c(.5, .5,.5)
y1 <- rdirichlet(50, omega1)

omega2 <- c(1., 1., 1.)
y2 <- rdirichlet(50, omega2)

omega3 <- c(1.5, 1.5, 1.5)
y3 <- rdirichlet(50, omega3)

y <- rbind(y1, y2, y3)

remove(list = c('omega1', 'y1', 'omega2', 'y2', 'omega3', 'y3'))

setwd('C:/Users/edoar/Desktop/Tesi/Code/BNP-change-point-detection-compositional')

n_iter <- 1000
burnin <- 500
q <- .5
trunc <- 3
iter_omega <- 100
burnin_omega <- 50
sigma_proposal_omega <- 1
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
                 iter_omega, burnin_omega, sigma_proposal_omega,
                 alpha_omega, beta_omega, 
                 alpha_sigma, beta_sigma,
                 alpha_propose_sigma, beta_propose_sigma,
                 alpha_theta, beta_theta)
