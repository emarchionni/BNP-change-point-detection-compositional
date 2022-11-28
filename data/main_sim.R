library('extraDistr')
n_data <- 100

omega <- c(0.7, 0.5, 1)

omega_norm <- sum(omega)

d <- length(omega)

trunc <- 100
t <- 1


y <- array(0, c(n_data, d))

y_old <- rdirichlet(1, omega)

y[1, ] <- y_old


for(datum in 2:n_data){
  
  m <- simulate_death_process(trunc, t, omega_norm)
  print(m)
  if(m == 0){
    
    l_vect <- rep(0, d) 
    
  } else {
    
    l_vect <- rmnom(1, m, y_old)
    
  }
    
  omega_sim <- omega + l_vect
  
  y_new <- rdirichlet(1, omega_sim)
  
  y[datum, ] <- y_new
  
  y_old <- y_new
  
  
}