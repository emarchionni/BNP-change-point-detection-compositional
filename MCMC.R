#'@param q: probability to perform a split
#'@param y: data, time X features
#'@param trunc: truncation value transition densities
#'@param step_omega: number of MH steps or number of samples in MC integration
#'@param MC_integration: if true, MCMC split&merge is performed using the integrated marginal likelihood, 
#'otherwise we perform an MH step for the value to estimate parameters of the the integrated likelihood,
#'note that in this latter case sigma_proposal_omega has no effect
#'@param save_data: if true, the output list will contain the data


# library('extraDistr')
# library('copula')
# library('Rcpp')

Rcpp::sourceCpp('transition_densities_omega_update_cpp.cpp')
source('MCMC_split_merge_shuffle.R')
source('MCMC_alpha.R')
source('MH_omega.R')
source('likelihood.R')
source('MC_likelihood.R')
source('eppf.R')
source('support_functions.R')
source('update_parameters/MCMC_sigma.R')
source('update_parameters/sample_theta.R')




MCMC <- function(niter, burnin, y, q,
                 trunc,
                 iter_omega, sigma_proposal_omega,
                 alpha_omega, beta_omega, 
                 alpha_sigma, beta_sigma,
                 alpha_propose_sigma, beta_propose_sigma,
                 alpha_theta, beta_theta,
                 MC_integration = F,
                 save_data = F){
  
  

  begin_time <- Sys.time()

  #### INITIALIZATION ####

  function_parameters <- as.list(environment(), all=TRUE)

  if(!save_data)
    function_parameters['y'] <- NULL

  # saving containers
  Acc_sigma <- array(0, niter)
  acc_sigma <- 0

  Theta <- array(0, niter)

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


  if(!MC_integration){
    # if MH step

    # initial omega
    omega <- array(rgamma(d * k, alpha_omega, beta_omega), c(k, d)) # cluster x component
    #omega <- array(1, c(k, d))

    # initial cluster likelihood
    likelihood <- full_log_integrated_likelihood(y, rho, omega, k, trunc)

  } else {
    # if MC integration

    # initial cluster likelihood
    likelihood <- MC_full_log_integrated_MARGINAL_likelihood(y, rho, trunc, n_clust,
                                                             iter_omega,
                                                             alpha_omega, beta_omega)

  }



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



    if(iter %% 50 == 0){
      print(iter)
      print(rho)
    }


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


      if(!MC_integration){
        # if MH step

        # update omega for new clusters
        omega_proposed <- MH_omega_split(iter_omega,
                                         y, j, d, rho_proposed,
                                         omega, sigma_proposal_omega,
                                         alpha_omega, beta_omega, trunc)

        # compute likelihood proposed
        likelihood_proposed <- full_log_integrated_likelihood_after_split(likelihood, y,
                                                                          rho_proposed, j, trunc,
                                                                          omega_proposed)


      } else {
        # if MH integration

        # MC integration likelihood
        likelihood_proposed <- MC_full_log_integrated_MARGINAL_likelihood_after_split(likelihood, y,
                                                                                      rho_proposed, j,
                                                                                      trunc,
                                                                                      iter_omega,
                                                                                      alpha_omega, beta_omega)

      }


      # compute log MH-alpha
      log_ratio <- MC_log_alpha_split(q, j,
                                      likelihood, eppf,
                                      likelihood_proposed, eppf_proposed,
                                      rho, rho_proposed)
      if(!is.double(log_ratio)|| is.na(exp(log_ratio) == 0)) browser() # TODO


      # MH step
      tot_partition_split <- tot_partition_split + 1


      if((log(runif(1)) <= min(0, log_ratio))){

        rho <- rho_proposed
        if(!MC_integration) omega <- omega_proposed
        likelihood <- likelihood_proposed
        eppf <- eppf_proposed
        acc_partition_split <- acc_partition_split + 1


      }





    } else {
      # merge

      # propose merge
      output_list <- merge(rho)

      rho_proposed <- output_list[[1]]
      j <- output_list[[2]]

      # compute eppf proposed
      eppf_proposed <- log_EPPF(rho_proposed, theta, sigma)


      if(!MC_integration){
        # if MH step

        # update omega for new clusters
        omega_proposed <- MH_omega_merge(iter_omega,
                                         y, j, rho_proposed,
                                         omega, sigma_proposal_omega,
                                         alpha_omega, beta_omega, trunc)


        # compute likelihood proposed
        likelihood_proposed <- full_log_integrated_likelihood_after_merge(likelihood, y,
                                                                          rho_proposed, j,
                                                                          trunc, omega_proposed)


      } else {
        # if MC integration

        # MC integration likelihood
        likelihood_proposed <- MC_full_log_integrated_MARGINAL_likelihood_after_merge(likelihood, y,
                                                                                      rho_proposed, j,
                                                                                      trunc,
                                                                                      iter_omega,
                                                                                      alpha_omega, beta_omega)

      }





      # compute log MH-alpha
      log_ratio <- MC_log_alpha_merge(q, j,
                                      likelihood, eppf,
                                      likelihood_proposed, eppf_proposed,
                                      rho, rho_proposed)
      if(!is.double(log_ratio)|| is.na(exp(log_ratio) == 0)) browser() #TODO


      # MH step
      tot_partition_merge <- tot_partition_merge + 1

      if((log(runif(1)) <= min(0, log_ratio))){

        rho <- rho_proposed
        if(!MC_integration) omega <- omega_proposed
        likelihood <- likelihood_proposed
        eppf <- eppf_proposed
        acc_partition_merge <- acc_partition_merge + 1

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


        if(!MC_integration){
          # if MH step

          # update omega for new clusters
          omega_proposed <- MH_omega_shuffle(iter_omega,
                                             y, j, d,
                                             rho, rho_proposed,
                                             omega, sigma_proposal_omega,
                                             alpha_omega, beta_omega, trunc)


          # compute likelihood proposed
          likelihood_proposed <- full_log_integrated_likelihood_after_shuffle(likelihood, y,
                                                                              rho, rho_proposed, j,
                                                                              trunc, omega_proposed)


        } else {
          # if MC integration

          # MC integration likelihood
          likelihood_proposed <- MC_full_log_integrated_MARGINAL_likelihood_after_shuffle(likelihood, y,
                                                                                          rho, rho_proposed, j,
                                                                                          trunc,
                                                                                          iter_omega,
                                                                                          alpha_omega, beta_omega)

        }




        # compute log MH-alpha


        log_ratio <- MC_log_alpha_shuffle(q, j,
                                          likelihood, eppf,
                                          likelihood_proposed, eppf_proposed,
                                          rho, rho_proposed)

        if(!is.double(log_ratio)|| is.na(exp(log_ratio) == 0)) browser() #TODO
        # MH step

        tot_partition_shuffle <- tot_partition_shuffle + 1

        if(log(runif(1)) <= min(0, log_ratio)){

          rho <- rho_proposed
          if(!MC_integration) omega <- omega_proposed
          likelihood <- likelihood_proposed
          eppf <- eppf_proposed
          acc_partition_shuffle <- acc_partition_shuffle + 1

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


  end_time <- Sys.time()

  return(list(call = function_parameters,
              partitions = Acc_partition_iter,
              theta = Theta,
              sigma = Acc_sigma,
              prop_acc_partition_split = (acc_partition_split / tot_partition_split),
              prop_acc_partition_merge = (acc_partition_merge / tot_partition_merge),
              prop_acc_partition_shuffle = (acc_partition_shuffle / tot_partition_shuffle),
              execution_time = (end_time - begin_time)))
}