paa <- function(list_coef, test_sim, n_training_sample, vars, r2_thres){

    Y_test = test_sim[["Y_list"]][[1]]
    Y_innov = test_sim[["innov_list"]][[1]]
    innov_maha <- colSums(t(Y_innov %*% solve(VAR_sigma)) * t(Y_innov))
    innov_maha <- innov_maha[2:length(innov_maha)]
    
    # Name variables
    var_names <- sprintf("V%d",seq(1:vars))
    lag_names <- sprintf("V%dlag",seq(1:vars))

    # Create an empty matrix to store the results of each model's predictive accuracy
    mat_pred_acc <- matrix(nrow = n_training_sample, ncol = length(r2_thres))
    r2_score_vec = c()

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
        pred_maha <- colSums(t(t(error_mat) %*% solve(VAR_sigma)) * error_mat)
        pred_maha <- pred_maha[2:length(pred_maha)]
        
        # Calculate predictive accuracy of an estimated model
        r2_score = cor(innov_maha, pred_maha)^2
        r2_score_vec = c(r2_score_vec, r2_score)

        for (i in 1:length(r2_thres)){
            thres = r2_thres[i]
            mat_pred_acc[p, i] <- r2_score > thres
        }
    }

    ### Step 4: Compute expected predictive accuracy ----
    PAP <- apply(mat_pred_acc, 2, function(x) sum(x) / n_training_sample)
    names(PAP) = paste0("PAP_",r2_thres)

    return(list(PAP=PAP, r2_score_vec=r2_score_vec))
}
