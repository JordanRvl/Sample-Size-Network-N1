pow <- function(list_pval, n_training_sample, vars, alpha){

    # Name variables
    var_names <- sprintf("V%d",seq(1:vars))
    lag_names <- sprintf("V%dlag",seq(1:vars))

    # Create dataframe to store power values of each parameters
    df_pow <- data.frame()

    for (var in 1:vars){
        # Get variable name
        name_ <- var_names[var]
        
        # Compute power
        df <- data.frame(t(colSums(list_pval[[var]] < alpha)/n_training_sample))

        # Rename df power
        names(df) <- c(paste0("int_", name_), paste0(name_, "_", lag_names))

        # Export
        df_pow <- bind_rows(df_pow, df)
    }

    # Keep only the non-missing value (1 per column)
    df_pow <- df_pow %>%
        summarise_all(~ ifelse(all(is.na(.)), NA, first(na.omit(.))))

    # Add "pow_" to the column names
    names(df_pow) <- paste0("pow_", names(df_pow))

    return(df_pow)
}
