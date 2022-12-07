### CALIBRATION MC integration ###

# 3-dim simplex



remove(list = ls())


library('foreach')
library('doParallel')

setwd('C:/Users/edoar/Desktop/Tesi/Code/BNP-change-point-detection-compositional')
#source('C:/Users/edoar/Desktop/Tesi/Code/BNP-change-point-detection-compositional/data/death_process.R')


path <- 'simulations/MC_calibration'

# grid of hyperparams to test
#hyperparams <- expand.grid(c(1), c(1))
hyperparams <- expand.grid(c(1, 2, 3), c(1, 20, 100, 250, 500))
colnames(hyperparams) <- c('trunc', 'iter_omega')
write.table(hyperparams, file = paste0(path, '/hyperparameters.txt'))
n_couples_param <- dim(hyperparams)[1]



# hyperparams to generate data
n_rep <- 5
n_change_points <- 3
n_data <- 50 # number of observation for each regime



omega_all <- rbind(c(0.2, 0.5, 1),
                   c(1, 0.5, 0.5),
                   c(0.3, 0.1, 0.2))
colnames(omega_all) <- c('1', '2', '3')
write.table(omega_all, file = paste0(path, '/omega.txt'))
omega_all/apply(omega_all, 1, sum)
trunc_sim <- 50
t <- 1
d <- dim(omega_all)[2]



# non-varying hyperparams of the algorithm
n_iter <- 4000
burnin <- 1000
q <- .5
sigma_proposal_omega <- NA
alpha_omega <- 2 
beta_omega <- 2
alpha_sigma <- 1
beta_sigma <- 1
alpha_propose_sigma <- 1 
beta_propose_sigma <- 1
alpha_theta <- 1
beta_theta <- 1


# register parallel backend
n_cores <- detectCores()
my_cluster <- makeCluster(n_cores-2)
registerDoParallel(cl = my_cluster)


# library('extraDistr')
# library('copula')
# library('Rcpp')


packages <- c('extraDistr',
              'Rcpp',
              'copula')

# output_list <- foreach(rep = 1:2,
#                        .inorder = FALSE,
#                        # .noexport = c("MCMC", )
#                        # .export = c("rdirichlet", "Stirling1", "rmnom"),
#                        .packages = packages,
#                        .verbose = T
# ) %dopar%
#   {
# 
#     Rcpp::sourceCpp('transition_densities_omega_update_cpp.cpp')
#   }


# parallelize over reps

output_list <- foreach(rep = 1:5, 
                       .inorder = FALSE,
                       .noexport = c("MCMC"),
                       # .export = c("rdirichlet", "Stirling1", "rmnom"),
                       .packages = packages,
                       .verbose = F
                       ) %dopar% 
  {
    # load functions
    setwd('C:/Users/edoar/Desktop/Tesi/Code/BNP-change-point-detection-compositional')
    source('C:/Users/edoar/Desktop/Tesi/Code/BNP-change-point-detection-compositional/data/death_process.R')
    source('MCMC.R')
    
    # log file
    sink(paste0(path, "/log.txt"), append=TRUE)
    
    
    # declare data structure
    y <- array(0, c(0, d))

    # GENERATE DATA FOR THE SINGLE REP

    for(point in 1:n_change_points){

      y_regime <- array(0, c(n_data, d))

      omega <- omega_all[point,]
      omega_norm <- sum(omega)

      y_old <- rdirichlet(1, omega)

      y_regime[1, ] <- y_old


      for(datum in 2:n_data){

        m <- simulate_death_process(trunc_sim, t, omega_norm)
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

    write.table(y, file = paste0(path,'/data_rep_', as.character(rep),'.txt'))

    list_curr_rep <- vector('list', n_couples_param)

    for(idx in 1:n_couples_param){

      print(paste('Replica', rep))
      print('Hyperparameters')
      print(paste('trunc:', hyperparams[idx, 1], 'nodes:', hyperparams[idx, 2]))

      trunc <- hyperparams[idx, 1]
      iter_omega <- hyperparams[idx, 2]

      #print(trunc)
      #print(iter_omega)

      list_curr_rep[[idx]] <- MCMC(n_iter, burnin, y, q,
                                   trunc,
                                   iter_omega, sigma_proposal_omega,
                                   alpha_omega, beta_omega,
                                   alpha_sigma, beta_sigma,
                                   alpha_propose_sigma, beta_propose_sigma,
                                   alpha_theta, beta_theta,
                                   TRUE, FALSE)


    }

    list_curr_rep

}


stopCluster(my_cluster)
