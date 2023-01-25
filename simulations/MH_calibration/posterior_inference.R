library(coda)
library(mcclust)
library(mcclust.ext)
library(foreach)
library(doParallel)

remove(list = ls())

setwd('C:/Users/edoar/Desktop/Tesi/Code/BNP-change-point-detection-compositional/simulations/MH_calibration')
load('output_MH_calibration.RData')
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


VI <- c(0, n_couples_param)
ESS_entropy <- c(0, n_couples_param)


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
