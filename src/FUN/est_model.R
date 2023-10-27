est_model <- function(name, n_training_sample, vars, Y_analysis){

    # Name variables
    var_names <- sprintf("V%d",seq(1:vars))
    lag_names <- sprintf("V%dlag",seq(1:vars))

    # List of dataframes to store estimates and pval
    list_coef <- list()
    list_pval <- list()
    for(var in 1:vars){
        list_coef[[var]] <- data.frame(matrix(nrow = n_training_sample, ncol = 1 + nrow(VAR_delta)))
        list_pval[[var]] <- data.frame(matrix(nrow = n_training_sample, ncol = 1 + nrow(VAR_delta)))
    }

    for (i in 1:n_training_sample){
         
        # Estimation exception for epskamp dataset
        if (name == "epskamp"){
            # Fit models
            fit_mod_V1 <- summary(lm(V1 ~ 1, data = Y_analysis[[i]]))
            fit_mod_V2 <- summary(lm(V2 ~ 1, data = Y_analysis[[i]]))
            fit_mod_V3 <- summary(lm(V3 ~ V7lag, data = Y_analysis[[i]]))
            fit_mod_V4 <- summary(lm(V4 ~ 1, data = Y_analysis[[i]]))
            fit_mod_V5 <- summary(lm(V5 ~ V5lag, data = Y_analysis[[i]]))
            fit_mod_V6 <- summary(lm(V6 ~ V6lag + V7lag, data = Y_analysis[[i]]))
            fit_mod_V7 <- summary(lm(V7 ~ V6lag, data = Y_analysis[[i]]))

            # Store coefficients      
            list_coef[[1]][i,] <- c(fit_mod_V1$coefficients[1,1], rep(0, 7)) 
            list_coef[[2]][i,] <- c(fit_mod_V2$coefficients[1,1], rep(0, 7))
            list_coef[[3]][i,] <- c(unname(fit_mod_V3$coefficients[1,1]), rep(0, 6), unname(fit_mod_V3$coefficients[2,1]))
            list_coef[[4]][i,] <- c(fit_mod_V3$coefficients[1,1], rep(0, 7))
            list_coef[[5]][i,] <- c(unname(fit_mod_V5$coefficients[1,1]), rep(0, 4), unname(fit_mod_V5$coefficients[2,1]), rep(0, 2))
            list_coef[[6]][i,] <- c(unname(fit_mod_V6$coefficients[1,1]), rep(0, 5), unname(fit_mod_V6$coefficients[2,1]), unname(fit_mod_V6$coefficients[3,1]))
            list_coef[[7]][i,] <- c(unname(fit_mod_V7$coefficients[1,1]), rep(0, 5), unname(fit_mod_V7$coefficients[2,1]), 0)
            
            # Store pvalues
            list_pval[[1]][i,] <- c(unname(fit_mod_V1$coefficients[1,4]), rep(1, 7))
            list_pval[[2]][i,] <- c(unname(fit_mod_V2$coefficients[1,4]), rep(1, 7)) 
            list_pval[[3]][i,] <- c(unname(fit_mod_V3$coefficients[1,4]), rep(1, 6), unname(fit_mod_V3$coefficients[2,4]))
            list_pval[[4]][i,] <- c(unname(fit_mod_V4$coefficients[1,4]), rep(1, 7))
            list_pval[[5]][i,] <- c(unname(fit_mod_V5$coefficients[1,4]), rep(1, 4), unname(fit_mod_V5$coefficients[2,4]), rep(1, 2))
            list_pval[[6]][i,] <- c(unname(fit_mod_V6$coefficients[1,4]), rep(1, 5), unname(fit_mod_V6$coefficients[2,4]), unname(fit_mod_V6$coefficients[3,4]))
            list_pval[[7]][i,] <- c(unname(fit_mod_V7$coefficients[1,4]), rep(1, 5), unname(fit_mod_V7$coefficients[2,4]), 0)

        } else {
            for (var in 1:vars){
                # Get variable name
                name_ <- var_names[var]

                # Create formula
                formula <- as.formula(paste(name_ ,paste(lag_names, collapse = " + "), sep = " ~ "))
                sum <- summary(lm(formula, data = Y_analysis[[i]])) 

                # Extract estimates
                list_coef[[var]][i,] <- unname(sum$coefficients[,1])

                # Extract pvalues
                list_pval[[var]][i,] <- unname(sum$coefficients[,4])
            }
        }
    }
    
    return(list("list_coef" = list_coef, "list_pval" = list_pval))
}
