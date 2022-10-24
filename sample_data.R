### sample data

library('extraDistr')
omega <- c(1, 1, 1)
y <- rdirichlet(50, omega)

n <- dim(y)[1]
d <- dim(y)[2]

rho <- c(as.integer(n/2))
rho <- c(rho, n-rho)




output <- split(rho)
rho_proposed <- output[[1]]
j <- output[[2]]
remove(output)

# y_0 <- y[1,]; y <- y[2,]


#split
omega_new <- array(0, c(3,d))
omega_new[1, ] <- c(1,1,1)
omega_new[2,] <- c(1.2,1.2,1.2)
omega_new[3,] <- c(1.5,1.5,1.5)
old_likelihood <- c(0,0)
omega <- omega_new


y_split <- split_data_partition(y, rho_proposed)

#sigma_0 omega 0.01
#MH_omega(100,500,y_split[[1]], omega[1,], 0.01, 1,1,16,10)

omega <- omega[1:2,]
MH_omega_split(100, 500, y, j, rho, rho_proposed, omega, 0.01, 1, 1, 5)
MH_omega_split(100, 200, y, j, rho, rho_proposed, omega, 0.01, 1, 1, 3)