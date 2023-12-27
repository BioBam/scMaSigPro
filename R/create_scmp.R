#' @title Create ScMaSigPro Object
#'
#' @description
#' `create_scmp()` initializes a ScMaSigPro object with the given counts,
#' cell level metadata, and other optional parameters.
#'
#' @param counts A matrix containing the raw expression counts.
#' @param cell_data A data frame containing the cell level metadata.
#' @param feature_data A data frame containing the feature level metadata.
#' @param use_as_bin A logical indicating to run MaSigPro analysis without
#' pseudotime binning. Defaults to FALSE.
#' @param bin_counts A matrix containing the binned counts. (Only if
#' use_as_bin == TRUE)
#' @param bin_cell_data A data frame containing the binned cell level metadata.
#' (Only if use_as_bin == TRUE)
#' @param ptime_col A character string representing the column name
#' for inferred Pseudotime values. (Default is "Pseudotime")
#' @param path_col A character string representing the column name for branching
#' path assignment. (Default is `path_prefix`)
#'
#' @return An object of class \code{\link{ScMaSigPro}}, with aligned pseudotime.
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#' @export

# Create ScMaSigPro
create_scmp <- function(counts,
                        cell_data,
                        feature_data,
                        bin_counts = NULL,
                        bin_cell_data = NULL,
                        ptime_col,
                        path_col,
                        use_as_bin = FALSE) {
  # Validation Checks
  assert_that(ncol(counts) == nrow(cell_data),
    msg = paste("Number of cells in raw-counts and cell-level-metadata are different.")
  )

  assert_that(all(colnames(counts) == rownames(cell_data)),
    msg = paste("Rownames of raw-counts and cell-level-metadata are different.")
  )

  if (!is.null(bin_counts) || !is.null(bin_cell_data)) {
    assert_that(nrow(bin_counts) == nrow(bin_cell_data),
      msg = paste("Number of cells in bin_counts and bin_cell_data are different.")
    )
    assert_that(all(rownames(bin_counts) == rownames(bin_cell_data)),
      msg = paste("Rownames of bin_counts and bin_cell_data are different.")
    )
  }

  # Create Single-Cell Experiment Object
  sparse_tmp <- SingleCellExperiment(
    list(counts = counts),
    colData = DataFrame(cell_data)
  )

  # Initate scMaSigPro
  scmpObj <- new("ScMaSigPro",
    Sparse = sparse_tmp,
    Dense = SingleCellExperiment(assays = list(bulk.counts = matrix(0, nrow = 0, ncol = 0)))
  )
  sparse_tmp <- NULL

  # Use as bin
  if (use_as_bin) {
    # Set bin size coulmn name
    cell_data[["scmp_bin_size"]] <- as.numeric(1)

    # Create sparse
    sparse_tmp <- SingleCellExperiment(
      list(bulk.counts = counts),
      colData = DataFrame(cell_data)
    )

    # Transfer Data
    scmpObj@Dense <- sparse_tmp
    sparse_tmp <- NULL

    # Update the slots
    scmpObj@Parameters@bin_ptime_col <- ptime_col
  }

  # Update the slots
  scmpObj@Parameters@ptime_col <- ptime_col
  scmpObj@Parameters@path_col <- path_col


  return(scmpObj)
}
