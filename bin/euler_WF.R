### EULER MARUYAMA METHOD FOR WF DIFFUSION ###

data_simulation_Wright_Fisher <- function(omega, n_points){
  
  # this fuction simulates equally spaced trajectories for a Wright-Fisher diffusion of parameter omega
  
  d <- length(omega)
  omega_norm <- sum(omega)
  
  y <- array(0, c(n_points, d))
  
  y_old <- rdirichlet(1, omega)
  y_new <- array(0, d)
  
  for (i in 2:n_points) { 
    # for each y to be simu,lated
    
    wiener <- rnorm(d - 1, 0, 1)
    
    for(j in 1:(d-1)){ 
      # for each component of y_new
      
      y_new[j] = y_old[j] + 0.5 * (omega[j] - omega_norm * y_old[j])
      
      for(k in 1:(d-1)){ 
        # for each component of y_old (for diffusion term)
        y_new[j] = y_new[j] + sqrt(y_old[j] * (1 - y_old[k])) * wiener[k]
        
      }
      
    }
    
    y_new[d] <- 1-sum(y_new[1:(d - 1)])
    
    y[i, ] <- y_new
    
    y_old <- y_new
    
  }
  
  return(y)
  
}
