
estimates_to_df = function(name, delta, phi, sigma){
    unroll_vec = function(mat, sig_mat=FALSE) {
        colnames(mat) = rownames(mat) = paste0("Y", 1:ncol(mat))
        values <- as.vector(mat) # Unroll into a named vector
        names <- paste(rep(rownames(mat), length(colnames(mat))), rep(colnames(mat), each = length(rownames(mat))), sep = "_")
        if (sig_mat) names <- paste0("sig_", names)
        setNames(values, names) # Create the named vector
    }

    name_data <- name
    names(name_data) <- "data"
    delta <- as.vector(delta)
    names(delta) <- paste0("int_Y", 1:length(delta))
    phi <- unroll_vec(phi)
    sigma <- unroll_vec(sigma, sig_mat=TRUE)

    data.frame(t(c(name_data, delta, phi, sigma)))
}
