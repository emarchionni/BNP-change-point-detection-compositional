n_params <- length(output_list)

total_time_rep <- 0

for(i in 1:n_params){
  total_time_rep <- total_time_rep + output_list[[i]][['1']][['execution_time']]
}



