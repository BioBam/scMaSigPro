#' Show ScMaSigPro Object Information
#'
#' This method displays basic information about the ScMaSigPro object when the object
#' is printed in the console. The method is automatically called when the user writes
#' the name of the object in the console.
#'
#' @param object An object of class \code{scMaSigProClass}.
#'
#' @importFrom S4Vectors coolcat
#'
#' @keywords internal
#' @export
.smsp_show <- function(object) {
  # Show Basic information
  cat("Class: ScMaSigPro\n")
  cat(paste0("nCells: ", ncol(object), "\n"))
  cat(paste0("nGenes: ", nrow(object), "\n"))
  cat("Continuum: ")
  cat(paste(range(object@covariate@time.series)[1],
    range(object@covariate@time.series)[2],
    sep = " to "
  ))
  cat("\n")

  # Calculate Dynamic Information
  if (length(object@pVector@adj.p.value) > 0) {
    sig.level <- object@parameters@p.vector.sig
    nSigs <- length(object@pVector@adj.p.value[object@pVector@adj.p.value <= sig.level])
    if (all(object@pVector@adj.p.value > sig.level)) {
      cat("Sig. Profiles (P-vector): None found")
    } else {
      cat(paste("Sig. Models (P-vector):", nSigs, sep = " "))
    }
  } else {
    cat("Sig. Models (P-vector): Not computed\n")
  }
}
