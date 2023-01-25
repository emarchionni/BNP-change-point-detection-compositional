vector_allocation <- function(partition){
  
  n_data <- sum(partition)
  n_clust <- length(partition)
  
  allocation_vector <- array(0, 0)
  
  clust_all <- 1
  
  for(clust in 1:n_clust){
    
    allocation_vector <- c(allocation_vector, 
                               rep(clust, partition[clust]))
    
    clust <- clust + 1
    
  }
  
  return(allocation_vector)
  
}


extract_allocation_matrix <- function(partition){
  
  if(length(partition) < 1)
    stop('Empty partition')
  
  n_sample <- length(partition)
  
  n_data <- sum(partition[[1]])
  
  allocation_matrix <- array(0, c(n_sample, n_data))
  
  for(i in 1:n_sample)
    allocation_matrix[i,] <- vector_allocation(partition[[i]])
  
  return(allocation_matrix)
  
}



extract_entropy <- function(partition){
  
  if(length(partition) < 1)
    stop('Empty partition')
  
  n_sample <- length(partition)
  
  n_data <- sum(partition[[1]])
  
  entropy <- c()
  
  for(i in 1:n_sample){
    
    curr_entropy <- 0
    
    curr_partition <- partition[[i]]
    
    k <- length(curr_partition)
    
    for(clust in 1:k)
      curr_entropy <- curr_entropy + curr_partition[clust]/n_data * log(curr_partition[clust]/n_data)
    
    
    entropy <- rbind(entropy, curr_entropy)
    
  }
  
  return(entropy)
}


# 
# 
# VI <- function(){
#   
# }