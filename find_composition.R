remove(list=ls())
K <- 3
N <- 30

n_vec <- array(0, c(0,K))


#total <- choose(N+K-1, N)



get_first_weak_composition <- function(n, k)
{
  composition <- array(0,c(1,K))
  if (n < k) {
    return(F)
  }
  
  for (i in 1:k-1) {
    composition[i] = 0
  }
  composition[k] = n
  #browser()
  return(composition)
}

get_next_weak_composition <- function(n, k, composition)
{
  if (composition[1] == n)    {
    return(F)
  }
  # there's an i with composition[i] > 0, and it is not 0. find the last one
    last = k
    
    while (composition[last] == 0) {
      #browser()
        last <- last-1
    }
    # turn    a b ...   y   z 0 0 ...   0
    #                      ^ last
    # into    a b ... (y+1) 0 0 0 ... (z-1)

    # be careful, there may be no 0's at the end
  
  z = composition[last]
  composition[last - 1] = composition[last - 1] + 1
  composition[last] = 0
  composition[k] = z - 1
  return(composition)
}




#### MAIN

composition <- get_first_weak_composition(N, K)


while(is.array(composition)){
  
  if(is.array(composition)){
    n_vec <- rbind(n_vec, composition)
  }
  
  composition <- get_next_weak_composition(N, K, composition)
  
  
}


