library(coda)
library(mcclust)
library(mcclust.ext)
library(foreach)
library(doParallel)

remove(list = ls())

setwd('C:/Users/edoar/Desktop/Tesi/Code/BNP-change-point-detection-compositional/simulations/MC_calibration')
load('output_MC_calibration.RData')
source("C:/Users/edoar/Desktop/Tesi/Code/BNP-change-point-detection-compositional/posterior inference/functions.R")

hyperparam <- read.table('hyperparameters.txt')
n_couples_param <- dim(hyperparam)[1]
n_rep <- length(output_list)

rho_true <- c(50, 50, 50)
rho_allocation_true <- vector_allocation(rho_true)

# register parallel backend
n_cores <- detectCores()
my_cluster <- makeCluster(n_cores-2)
registerDoParallel(cl = my_cluster)



### ALLOCATION PARTITIONS ###


rho_allocation <- foreach(rep = 1:5, 
                       .inorder = TRUE,
                       # .export = c(),
                       #.packages = packages,
                       .verbose = F
) %dopar% 
  {

    
    list_curr_allocation <- vector('list', n_couples_param)
    
    for(idx in 1:n_couples_param)
      list_curr_allocation[[idx]] <- extract_allocation_matrix(output_list[[rep]][[idx]][['partitions']])
  
    
    list_curr_allocation
    
  }



### ESTIMATED PARTITIONS ###

packages <- c('mcclust', 'mcclust.ext')

estimated_partitions <- foreach(rep = 1:5, 
                       .inorder = TRUE,
                       # .export = c("rdirichlet", "Stirling1", "rmnom"),
                       .packages = packages,
                       .verbose = F
) %dopar% 
  {
    
    list_estimated_allocation <- vector('list', n_couples_param)
    
    for(idx in 1:n_couples_param)
      list_estimated_allocation[[idx]] <- minVI(comp.psm(rho_allocation[[rep]][[idx]]), 
                                                rho_allocation[[rep]][[idx]], 
                                                method = "draws")
    
    list_estimated_allocation
      
    
  }



### VI ESTIMATED PARTITION - REAL PARTITION ###


packages <- c('mcclust')

metrics_VI <- foreach(rep = 1:5, 
                      .inorder = TRUE,
                      # .export = c(),
                      .packages = packages,
                      .verbose = F
) %dopar% 
  {
    
    metrics_curr <- vector('list', n_couples_param)
    
    for(idx in 1:n_couples_param)
      metrics_curr[[idx]][['VI']] <- vi.dist(rho_allocation_true, 
                                             estimated_partitions[[rep]][[idx]][['cl']])
    
    metrics_curr
    
    
  }



### ESS ENTROPY ESTIMATED PARTITION - REAL PARTITION ###


entropy <- foreach(rep = 1:5, 
                   .inorder = TRUE,
                   # .export = c(),
                   #.packages = packages,
                   .verbose = F
) %dopar% 
  {
    
    entropy_curr <- vector('list', n_couples_param)
    
    
    for(idx in 1:n_couples_param)
      entropy_curr[[idx]] <- extract_entropy(output_list[[rep]][[idx]][['partitions']])
    
    
    entropy_curr
    
    
  }


packages <- c('coda')

metrics_entropy <- foreach(rep = 1:5, 
                           .inorder = TRUE,
                           # .export = c(),
                           .packages = packages,
                           .verbose = F
) %dopar% 
  {
    
    metrics_entropy_curr <- vector('list', n_couples_param)
    
    
    for(idx in 1:n_couples_param)
      metrics_entropy_curr[[idx]][['ESS_entropy']] <- as.double(effectiveSize(entropy[[rep]][[idx]]))
    
    
    metrics_entropy_curr
    
    
  }


VI <- c()
ESS_entropy <- c()


for(idx in 1:n_couples_param){
  
  curr_mean_VI <- 0
  curr_mean_ESSentropy <- 0
  
  for(rep in 1:n_rep){
    
    curr_mean_VI <- curr_mean_VI + metrics_VI[[rep]][[idx]][['VI']]
    curr_mean_ESSentropy <- curr_mean_ESSentropy + metrics_entropy[[rep]][[idx]][['ESS_entropy']]
    
  }
  
  curr_mean_VI <- curr_mean_VI / n_rep
  curr_mean_ESSentropy <- curr_mean_ESSentropy / n_rep
  
  VI[idx] <- curr_mean_VI
  ESS_entropy[idx] <- curr_mean_ESSentropy
  
}

write.table(VI, file = 'posterior_inference/VI.txt')
write.table(ESS_entropy, file = 'posterior_inference/ESS_entropy.txt')


stopCluster(my_cluster)


### SAMPLING TIME ###


sampling_time <- c()


for(idx in 1:n_couples_param){
  
  curr_mean <- 0
  
  for (rep in 1:5){
    units(output_list[[rep]][[idx]][['execution_time']]) <- 'secs'
    curr_mean <- curr_mean + output_list[[rep]][[idx]][['execution_time']]
  }
    
    
  
  curr_mean <- curr_mean/n_rep
  sampling_time[idx] <- curr_mean

  
}

write.table(sampling_time, file = 'posterior_inference/sampling_time.txt')

#### PLOTS ####


library(ggplot2)

setwd('C:/Users/edoar/Desktop/Tesi/Code/BNP-change-point-detection-compositional/simulations/MC_calibration')


hyperparam <- read.table('hyperparameters.txt')
n_couples_param <- dim(hyperparam)[1]

rho_true <- c(50, 50, 50)

entropy <- read.table('posterior_inference/ESS_entropy.txt')
VI <- read.table('posterior_inference/VI.txt')
sampling_time <- read.table('posterior_inference/sampling_time.txt')

list_VI <- vector('list', 3)
list_entropy <- vector('list', 3)
list_sampling_time <- vector('list', 3)

for(trunc in 1:3){
  
  for(i in 0:4){
   
    list_VI[[trunc]][[i+1]] <- VI[trunc + i * 3, 1]
    list_entropy[[trunc]][[i+1]] <- entropy[trunc + i * 3, 1]
    list_sampling_time[[trunc]][[i+1]] <- sampling_time[trunc + i * 3, 1]
  }
  
}

omega_iter <- c(1,20,100,250,500)


### VI ###

# level of truncation 1

trunc <- 1

VI <- array(0, 5)
time <- array(0, 5)

for(i in 1:5){
  
  VI[i] <- list_VI[[trunc]][[i]]
  time[i] <- list_sampling_time[[trunc]][[i]]
  
  
}

curr_VI <- data.frame(omega_iter, VI, time)


ggplot(curr_VI, aes(x = omega_iter, y = VI)) + 
  geom_point(colour = '#FF9933', size = 2.5) +
  labs(x = 'Number of ierations',
       y = 'Variation of information',
       subtitle = 'Level of truncation: 1') + 
  theme_minimal()

# ggplot(curr_VI, aes(x = omega_iter, y = VI, colour=time)) + 
#   geom_point(size = 2.5) +
#   labs(x = 'Number of ierations',
#        y = 'Variation of information',
#        subtitle = 'Level of truncation: 1') + 
#   theme_minimal()



# level of truncation 2

trunc <- 2

VI <- array(0, 5)

for(i in 1:5){
  
  VI[i] <- list_VI[[trunc]][[i]]
  
}

curr_VI <- data.frame(omega_iter,VI)

ggplot(curr_VI, aes(x = omega_iter, y = VI)) + 
  geom_point(colour = 'orange', size = 2) +
  labs(x = 'Number of ierations',
       y = 'VI',
       subtitle = 'Level of truncation: 2') + 
  theme_minimal()


# level of truncation 3

trunc <- 3

VI <- array(0, 5)

for(i in 1:5){
  
  VI[i] <- list_VI[[trunc]][[i]]
  
}

curr_VI <- data.frame(omega_iter,VI)

ggplot(curr_VI, aes(x = omega_iter, y = VI)) + 
  geom_point(colour = 'orange', size = 2) +
  labs(x = 'Number of ierations',
       y = 'VI',
       subtitle = 'Level of truncation: 3') + 
  theme_minimal()


### ESS ###
