#' @importFrom S4Vectors coolcat
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

#' @export
setMethod("show", "scMaSigProClass", .smsp_show)
