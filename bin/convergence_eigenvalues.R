remove(list = ls())

nn <- seq(0,50,1)
omega <- c(0.05, 0.5, 1, 5, 10, 30)
#exp_lambda : (i,j) i: value of |n| in the vector nn j: omega considered
exp_lambda <- array(0, dim = c(length(nn),length(omega)))

for (i in 1:length(omega)) {
  exp_lambda[,i] <- exp(-0.5 * nn * (nn - rep(1,length(nn)) + omega[i]))
}
