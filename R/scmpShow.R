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
.scmp_show <- function(object) {
  # Show Basic information
  cat("Class: ScMaSigPro\n")
  cat(paste0("nCells: ", ncol(object@sce), "\n"))
  cat(paste0("nFeatures: ", nrow(object@sce), "\n"))
  cat("Pseudotime Range:", paste(range(colData(object@sce)[[object@addParams@pseudotime_colname]])))

  # Calculate the Compression
  compressed.cell.metadata <- object@compress.sce@colData %>% as.data.frame()
  if (length(compressed.cell.metadata) > 0) {
    cat(paste("\nPaths:", paste(levels(as.factor(compressed.cell.metadata[[object@addParams@path_colname]])), collapse = ", ")))
    cat(paste0(
      "\nBinned Pseudotime: ", paste(range(compressed.cell.metadata[[object@addParams@bin_pseudotime_colname]]), collapse = "-"), "(Range), ",
      round(mean(compressed.cell.metadata[[object@addParams@bin_pseudotime_colname]]), 2), "(Mean), "
    ))

    # Extract info
    per_path_num_bin <- extract_info(compressed.cell.metadata, return_type = "num_bins", bin_size_col = object@addParams@bin_size_colname, object@addParams@path_colname)
    per_path_bin_size <- round(extract_info(compressed.cell.metadata, return_type = "avg_bin_size", bin_size_col = object@addParams@bin_size_colname, object@addParams@path_colname))

    # Paste
    cat("\nNumber of bins->", paste(names(per_path_num_bin), per_path_num_bin, sep = ": "))
    cat("\nAverage bin Size->", paste(names(per_path_bin_size), per_path_bin_size, sep = ": "))
  }

  # Calculate Dynamic Information
  if (length(object@scPVector@p.adjusted) > 0) {
    sig.level <- object@addParams@Q
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

# helper to extract the lineage info
extract_info <- function(data, return_type = "avg_bin_size", bin_size_col, path_col) {
  if (return_type == "avg_bin_size") {
    avg_sizes <- tapply(data[[bin_size_col]], data[[path_col]], mean)
    return(avg_sizes)
  } else if (return_type == "num_bins") {
    bin_counts <- table(data[[path_col]])
    return(bin_counts)
  } else {
    stop("Invalid return_type. Choose between 'avg_bin_size' and 'num_bins'.")
  }
}
