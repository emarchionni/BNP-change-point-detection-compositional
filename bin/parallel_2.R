library('foreach')
library('doParallel')

n_cores <- detectCores()
my_cl <- makeCluster(n_cores-2)
registerDoParallel(cl = my_cl)




a <- foreach(i = 1:10, .packages = 'extraDistr') %dopar% {
  rdirichlet(i, c(0.2,0.3,1))
}
  

stopCluster(my_cl)
