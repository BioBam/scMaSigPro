showTscores <- function(scMaSigPro.obj, view = T) {
  # Extract Coeffcients
  t_value <- sapply(scMaSigPro.obj@tFit@model.attributes, simplify = T, function(x) {
    return(x[["t_value"]])
  })

  # Transpose
  t_value <- t(t_value)

  # If viewing is requested
  if (view == T) {
    View(t_value)
  }

  # Return
  return(t_value)
}
