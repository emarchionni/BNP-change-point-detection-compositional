
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

setwd('C:\\Users\\edoar\\Desktop\\Tesi\\Code\\BNP-change-point-detection-compositional\\testing\\transition densities')
library('Rcpp')
sourceCpp('transition_densities_cluster_cpp_testing.cpp')
source('transition_densities.R')
source('support_functions.R')

trunc <- 3
omega <- c(0.8,0.8,0.8)
log_transition_densities(y[0,], y[1,], omega, 3)
log_transition_densities_cpp(y[0,], y[1,], omega, 3)
system.time(log_transition_densities(y[0,], y[1,], omega, 3))
system.time(log_transition_densities_cpp(y[0,], y[1,], omega, 3))


source('likelihood.R')
n_clust <- dim(y)[1]
system.time(log_integrated_likelihood_cluster(y, omega, n_clust, trunc))
system.time(log_integrated_likelihood_cluster_cpp(y, omega, n_clust, trunc))
