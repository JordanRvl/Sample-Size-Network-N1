sim_var_data <- function(tm, n_training_sample, VAR_delta, VAR_phi, VAR_sigma){
    # Set the number of burning observations to be included
    N_burning <- 1000
    N_sim <- N_burning + tm

    # For one run:
    # Simulate innovations
    sim_innovation_list <- list()
    innovation_list <- list()

    for (i in 1:n_training_sample){
        sim_innovation <- mvrnorm(n = N_sim,
                                    mu = rep(0, nrow(VAR_delta)),
                                    Sigma = VAR_sigma)
        sim_innovation_list[[i]] <- sim_innovation
    }

    # Set the starting value of the time series as the expected value + innovation 
    # at the first timepoint
    VAR_expected_value <- solve(diag(nrow(VAR_delta)) - VAR_phi) %*% VAR_delta

    Y_start_list <- list()

    for (i in 1:n_training_sample){
        Y_start <- VAR_expected_value + sim_innovation_list[[i]][1,]
        Y_start_list[[i]] <- Y_start
    }


    # Simulate the entire time series
    Y_complete_list <- list()

    for (i in 1:n_training_sample){
        Y_complete <- matrix(0, nrow = N_sim, ncol = nrow(VAR_delta))
        Y_complete[1,] <- Y_start_list[[i]]
        for (t in 2: N_sim){
            Y_complete[t,] <- VAR_delta + 
            VAR_phi%*%Y_complete[t-1,] +  
            sim_innovation_list[[i]][t,]
        }
        Y_complete_list[[i]] <- Y_complete
    }

    # Exclude burning observations
    Y_list <- list()
    for (i in 1:n_training_sample){
        Y <- Y_complete_list[[i]][-(1:N_burning),]
        Y_list[[i]] <- Y

        innovation <- sim_innovation_list[[i]][-(1:N_burning),]
        innovation_list[[i]] <- innovation
    }

    # Create lagged variables for analysis and rename variables
    Y_analysis_list <- list()

    for (i in 1:n_training_sample){
        Y_analysis <- as.data.frame(cbind(Y_list[[i]], dplyr::lag(Y_list[[i]])))
        colnames(Y_analysis) <- c(paste0("V", 1:vars), sprintf("V%dlag",seq(1:vars)))
        Y_analysis_list[[i]] <- Y_analysis
    }
    
    return(list(Y_list=Y_analysis_list, innov_list=innovation_list))
}
