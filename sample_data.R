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

y_0 <- y[1,]; y <- y[2,]
