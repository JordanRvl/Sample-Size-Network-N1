paa <- function(list_coef, Y_test, n_training_sample, vars, maha_thres1, maha_thres2){
    
    # Name variables
    var_names <- sprintf("V%d",seq(1:vars))
    lag_names <- sprintf("V%dlag",seq(1:vars))

    # Create an empty matrix to store the results of each model's predictive accuracy
    mat_pred_acc <- matrix(nrow = n_training_sample, ncol = 1)

    # Transform Y_test
    Y_test_out <- Y_test[, var_names]
    Y_test_pred <- cbind(data.frame("int" = rep(1, nrow(Y_test))), Y_test[, lag_names])

    for (p in 1:n_training_sample){
        
        # Calculate the prediction error for each observation in the test set 
        coefs <- unlist(lapply(1:vars, function(x) list_coef[[x]][p, ]))
        VAR_mat_est <- matrix(coefs, byrow = T, nrow = vars)
        error_mat <- VAR_mat_est %*% t(Y_test_pred) - t(Y_test_out)
        
        # Calculate squared Mahalanobis distance from the error matrix
        # The slow way: calculate the entire matrix where values on diagonal are the target
        #Maha_vec <- diag(t(error_mat) %*% solve(VAR_sigma) %*% error_mat)
        # The fast way: only calculate the diagonal values in the matrix multiplication
        Maha_vec <- colSums(t(t(error_mat) %*% solve(VAR_sigma)) * error_mat)
        
        # Calculate predictive accuracy of an estimated model
        mat_pred_acc[p, 1] <- sum(Maha_vec < qchisq(maha_thres1, df = vars), na.rm=TRUE)/(nrow(Y_test))
    }

    # Compare it to the pre-defined performance threshold (use .94 here)
    N_good_model <- sum(mat_pred_acc >= maha_thres2)

    ### Step 4: Compute expected predictive accuracy ----
    PAP <- N_good_model/n_training_sample

    return(PAP)
}
