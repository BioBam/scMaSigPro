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
  cat(paste0("nCells: ", ncol(object@sce), "\n"))
  cat(paste0("nGenes: ", nrow(object@sce), "\n"))
  cat("Continuum: Toti ")
  cat("\n")

  # Calculate Dynamic Information
  if (length(object@scPVector@p.adjusted) > 0) {
    sig.level <- object@scPVector@Q
    nSigs <- length(object@scPVector@p.adjusted[object@scPVector@p.adjusted <= sig.level])
    if (all(object@scPVector@p.adjusted > sig.level)) {
      cat("Sig. Profiles (P-vector): None found")
    } else {
      cat(paste("Sig. Models (sc.p.vector):", nSigs, sep = " "))
    }
  } else {
    cat("Sig. Models (P-vector): Not computed\n")
  }
}
