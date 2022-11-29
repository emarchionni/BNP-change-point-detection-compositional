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

n_iter <- 1500
burnin <- 0
q <- .8
trunc <- 3
iter_omega <- 1000
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
               trunc,
               iter_omega, sigma_proposal_omega,
               alpha_omega, beta_omega, 
               alpha_sigma, beta_sigma,
               alpha_propose_sigma, beta_propose_sigma,
               alpha_theta, beta_theta, 
               TRUE)




#### RUN ####

remove(list = ls())


library('extraDistr')
library('MASS')
library('copula')
library('Compositional')

# data simulation for O-U process

gamma_sim_1 = 0.5
gamma_sim_2 = 0.5
gamma_sim_3 = 0.5

mu_1 <- c(-0.2,-0.4)
mu_2 <- c(1,1)
mu_3 <- c(-2,2)


sigma_1 = sigma_2 = sigma_3 = matrix(0, nrow = 2, ncol = 2)


diag(sigma_1) = c(0.5,0.2)

sigma_1[1,2]  = 0.2
sigma_1[2,1]  = 0.2



diag(sigma_2) <- c(1.5,2)

diag(sigma_3) <- c(2.5,2.5)
sigma_3[2,1]  = 0.4
sigma_3[1,2]  = 0.4

data_scenario_1 <- as.data.frame(matrix(nrow = 150, ncol = 2))

data_scenario_1[1,] = mu_1

for(i in 2:50){ 
  
  data_scenario_1[i,] = gamma_sim_1*data_scenario_1[i-1,] + (1-gamma_sim_1)*mu_1 + mvrnorm(n = 1,mu = mu_1, Sigma =  sigma_1)
  
}


data_scenario_1[51,] = mu_2

for(i in 52:100){ 
  
  data_scenario_1[i,] = gamma_sim_2*data_scenario_1[i-1,] + (1-gamma_sim_2)*mu_2 + mvrnorm(n = 1,  mu = mu_2, Sigma =  sigma_2)
  
}

data_scenario_1[100,] = mu_3

for(i in 101:150){ 
  
  data_scenario_1[i,] = gamma_sim_3*data_scenario_1[i-1,] + (1-gamma_sim_3)*mu_3 + mvrnorm(n = 1, mu = mu_3, Sigma =  sigma_3)
  
}

data <- as.matrix(data_scenario_1)

y <- as.matrix(alrinv(data))

dev.off()
plot(y[1:150,1],y[1:150,2])
plot(y[1:150,2],y[1:150,3])
plot(y[1:150,1],y[1:150,3])



setwd('C:/Users/edoar/Desktop/Tesi/Code/BNP-change-point-detection-compositional')

n_iter <- 30000
burnin <- 0
q <- .8
trunc <- 3
iter_omega <- 150
#burnin_omega <- 50
sigma_proposal_omega <- 3
alpha_omega <- .1
beta_omega <- .1
alpha_sigma <- 1
beta_sigma <- 1
alpha_propose_sigma <- 1 
beta_propose_sigma <- 1
alpha_theta <- 1
beta_theta <- 1


source('MCMC.R')

output <- MCMC(n_iter, burnin, y, q,
               trunc,
               iter_omega, sigma_proposal_omega,
               alpha_omega, beta_omega, 
               alpha_sigma, beta_sigma,
               alpha_propose_sigma, beta_propose_sigma,
               alpha_theta, beta_theta, 
               F)

save(output, file = "simulation_17_11_30Kiter_2.RData")
