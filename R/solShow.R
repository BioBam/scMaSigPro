showSol <- function(scMaSigPro.obj, view = T) {
  # Extract Coeffcients
  sol <- sapply(scMaSigPro.obj@tFit@model.attributes, simplify = T, function(x) {
    return(x[["sol"]])
  })

  # Transpose
  sol <- t(sol)

  # If viewing is requested
  if (view == T) {
    View(sol)
  }

  # Return
  return(sol)
}
