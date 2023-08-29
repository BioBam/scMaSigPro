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
  cat(paste0("nFeatures: ", nrow(object@sce), "\n"))
  cat("Continuum:")

  # Calculate the Compression
  compressed.cell.metadata <- as.data.frame(colData(object@compress.sce))
  if (length(compressed.cell.metadata) > 0) {
    cat(paste("\nPaths:", paste(levels(as.factor(compressed.cell.metadata$path)), collapse = ", ")))
    cat(paste0(
      "\nBinned Pseudotime: ", paste(range(compressed.cell.metadata$binnedTime), collapse = "-"), "(Range), ",
      mean(compressed.cell.metadata$bin.size), "(Mean), ",
      median(compressed.cell.metadata$bin.size), "(Median)"
    ))
  }

  # Calculate Dynamic Information
  if (length(object@scPVector@p.adjusted) > 0) {
    sig.level <- object@scPVector@Q
    nSigs <- length(object@scPVector@p.adjusted[object@scPVector@p.adjusted <= sig.level])
    if (all(object@scPVector@p.adjusted > sig.level)) {
      cat("\nSig. Profiles (P-vector): None found")
    } else {
      cat(paste("\nSig. Models (sc.p.vector):", nSigs, sep = " "))
    }
  }

  # Influential Genes if any
  if (ncol(object@scTFit@influ.info) > 0) {
    cat(paste("\nNo. of Influential Features:", ncol(object@scTFit@influ.info)))
  }
}
