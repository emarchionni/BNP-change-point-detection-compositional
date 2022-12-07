n.cores <- parallel::detectCores() - 1

my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")

doParallel::registerDoParallel(cl = my.cluster)

x <- foreach(
  i = 1:10, 
  .combine = 'c'
) %dopar% {
  sqrt(i)
}

parallel::stopCluster(cl = my.cluster)


registerDoSEQ() # back to sequential

