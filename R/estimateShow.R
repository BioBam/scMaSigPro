showEstimates <- function(scMaSigPro.obj, view = T){

    # Extract Coeffcients
    beta_estimates <- sapply(scMaSigPro.obj@tFit@model.attributes, simplify = T, function(x) {
        return(x[["beta_estimates"]])
        })

    # Transpose
    beta_estimates <- t(beta_estimates)

    # If viewing is requested
    if (view ==T){
        View(beta_estimates)
    }

    # Return
    return(beta_estimates)
}
