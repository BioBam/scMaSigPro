#' Align pseudotime from two different paths
#'
#' @description
#' `align.pseudotime()` is an internal function that aligns the pseudotime from
#' two different trajectory paths. It first identifies the path which has a longer
#' pseudotime range, and then it uses the range of the shorter path to rescale
#' the range of the longer one. It uses a simple rescaling mechanism from
#' `scales::rescale`. In future we expect to support strategies like dynamic time
#' warping from dtw package.
#'
#' @importFrom scales rescale
#' 
#' @param pseudotime_col Chacarcter string with the column name in `cell.metadata` storing
#' Pseudotime values. It is generated using `colData` from the \pkg{SingleCellExperiment}
#' package. (Default is "Pseudotime")
#' @param path_col Chacarcter string with the column name in `cell.metadata` storing information
#' for Path. It is generated using `colData` from the \pkg{SingleCellExperiment} package.
#' (Default is `path_prefix`)
#' @param method Currently only `scales::rescale` is supported. (Default is "rescale")
#' @param verbose Print detailed output in the console. (Default is TRUE)
#'
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#' @keywords internal

align.pseudotime <- function(scmpObj, pseudotime_col, path_col, method = "rescale", verbose = TRUE) {
  # Extract Cell metadata
  cell.metadata <- scmpObj@sce@colData %>% as.data.frame()

  # Extract Time and Path info with cell
  cell.metadata.sub <- cell.metadata[, c(pseudotime_col, path_col), drop = FALSE]
  cell.metadata.sub$cell <- rownames(cell.metadata.sub)

  # Get paths
  path.vec <- unique(cell.metadata.sub[[path_col]])

  # Get paths
  path1_time <- cell.metadata.sub[cell.metadata.sub[[path_col]] == path.vec[1], pseudotime_col]
  names(path1_time) <- cell.metadata.sub[cell.metadata.sub[[path_col]] == path.vec[1], "cell"]
  path2_time <- cell.metadata.sub[cell.metadata.sub[[path_col]] == path.vec[2], pseudotime_col]
  names(path2_time) <- cell.metadata.sub[cell.metadata.sub[[path_col]] == path.vec[2], "cell"]

  # Get list
  pTimeVectors <- select_longer_vector(
    vector1 = path1_time, vector1_label = path.vec[1],
    vector2 = path2_time, vector2_label = path.vec[2]
  )

  if (length(pTimeVectors) == 4) {
    if (verbose) {
      message(paste(pTimeVectors$long_vec_label, "is greater than", pTimeVectors$short_vec_label))
      message(paste("Adjusting", pTimeVectors$long_vec_label, "from", paste(pTimeVectors$short_vec_label, collapse = "-")))
    }

    # Check for requested method
    if (method == "rescale") {
      # Update
      pTimeVectors$long_vec <- rescale(pTimeVectors$long_vec, to = c(min(pTimeVectors$short_vec), max(pTimeVectors$short_vec)))

      if (verbose) {
        message(paste("New range for the pseudotime of", pTimeVectors$long_vec_label, "is now rescaled tp", paste(range(pTimeVectors$short_vec), collapse = "-")))
      }

      short_tmp <- data.frame(
        time = pTimeVectors$short_vec,
        cell = names(pTimeVectors$short_vec)
      )
      colnames(short_tmp) <- c(paste(pseudotime_col, "rescaled", sep = "_"), "cell")
      long_tmp <- data.frame(
        time = pTimeVectors$long_vec,
        cell = names(pTimeVectors$long_vec)
      )
      colnames(long_tmp) <- c(paste(pseudotime_col, "rescaled", sep = "_"), "cell")
      new_time <- rbind(short_tmp, long_tmp)

      # Merge
      cell.metadata <- merge(cell.metadata, new_time, by = "cell")
      rownames(cell.metadata) <- cell.metadata[["cell"]]

      # Drop cell
      cell.metadata <- cell.metadata[, colnames(cell.metadata) != "cell", drop = FALSE]

      # Add
      scmpObj@sce@colData <- DataFrame(cell.metadata)

      # Update Pseudotime
      scmpObj@addParams@pseudotime_colname <- paste(pseudotime_col, "rescaled", sep = "_")
    } else {
      if (verbose) {
        message("Currently only 'scales::rescale' is supported")
      }
      return(scmpObj)
    }
    return(scmpObj)
  } else {
    if (verbose) {
      message("Both paths are comparable, pseudotime is not rescaled")
    }
    return(scmpObj)
  }
}
