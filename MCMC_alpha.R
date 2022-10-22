#### ALPHA ####



MC_alpha <- function(trunc, step_omega, y, omega, j, rho_proposed, rho){
  
  k_proposed <- length(rho_proposed)
  k_previous <- length(rho)
  
  y_proposed <- split_data_partition(y, rho_proposed)
  y_previous <- split_data_partition(y, rho)
  
  #we pass old_omega, we compute new omega  
    

  
}


