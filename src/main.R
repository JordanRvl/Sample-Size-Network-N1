# Package loading
library(dplyr)
library(tidyr)
library(qgraph)
library(MASS)
# library(foreign)
library(graphicalVAR)
library(tictoc)


# Renv management
# renv::init()
# renv::snapshot()
# renv::restore()

# Import custom functions in src/FUN folder
for (f in list.files("src/FUN")){
    source(paste0("src/FUN/", f))
}

# Run parameters estimation
# source("src/estimate_parameters.R")

# Load parameters
load("data/list_param.rda")




########################
######## Parameter simulation
########################

# Simulation parameters
sim_timepoints <- c(seq(25,250,25), 63, 86, 47, 300, 353, 400, 500) # Ttraining
sim_timepoints <- sim_timepoints[order(sim_timepoints)]
n_training_sample <- 5000 # TODO: change to 5000

# PAA parameters
test_timepoints <- 100000
r2_thres <- c(.9, .92, .95)

# Power parameters
alpha <- .05
power_threshold <- .8

# Seed: obtain using runif(1) * 10000000
set.seed(1921325)



########################
######## Simulation
########################

# Loop over each parameter of the simulation
for (name in names(list_param)){
    print(paste0("----- Parameters number:", name, " -----"))

    # Dataframe to store monte carlo simulation results
    df_mc <- data.frame()

    # Import parameters
    vars <-  list_param[[name]][["vars"]]
    VAR_delta <- list_param[[name]][["delta"]]
    VAR_phi <- list_param[[name]][["phi"]]
    VAR_sigma <-list_param[[name]][["sigma"]]

    # Loop over each T pre-defined above
    for (i in 1:length(sim_timepoints)){
        start <- Sys.time() 

        # Select timepoints training
        tm <- sim_timepoints[i]
        print(paste0("T = ", tm))
        tic()
        
        # Simulate Test set
        test_sim <- sim_var_data(test_timepoints, 1, VAR_delta, VAR_phi, VAR_sigma)

        # Simulate datasets
        Y_analysis <- sim_var_data(tm, n_training_sample, VAR_delta, VAR_phi, VAR_sigma)[["Y_list"]]

        # Estimate models
        list_model <- est_model(name, n_training_sample, vars, Y_analysis)

        # Estimate power 
        res_pow <- pow(list_pval=list_model$list_pval, n_training_sample, vars, alpha)

        # Estimate PAA
        res_paa <- paa(list_coef=list_model$list_coef, test_sim, n_training_sample, vars, r2_thres)

        # Get simulation time in seconds
        time <- round(as.numeric(Sys.time() - start, units="secs"), 2)
        print(paste0("Time = ", time))
        toc()

        # Store results in a dataframe
        df <- bind_cols(data.frame(name=name, vars=vars, timepoint=tm, sim_time = time),
                        data.frame(t(res_paa[["PAP"]])),
                        res_pow)
        df_mc <- bind_rows(df_mc, df)
    }

    # Save results
    if(!dir.exists("data/data_sim")) dir.create("data/data_sim")
    write.csv(df_mc, paste0("data/data_sim/df_results_", name,".csv"), row.names=FALSE)
}

