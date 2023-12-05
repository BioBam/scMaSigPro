#' @title Create scMaSigProClass Object
#'
#' @description
#' `create.scmp()` initializes a scMaSigProClass object with the given counts,
#' cell level metadata, and other optional parameters.
#'
#' @param counts A matrix containing the raw expression counts.
#' @param cell_data A data frame containing the cell metadata.
#' @param feature_data A data frame containing the feature level metadata.
#' @param bin_counts A matrix containing the binned counts.
#' @param bin_cell_data A data frame containing the binned cell level metadata.
#' @param pseudotime_colname A character string specifying the column name for Pseudotime in the cell level metadata.
#' @param path_colname A character string specifying the column name for the path in the cell level metadata.
#' @param use_as_bin A logical indicating whether to use the raw counts and cell level data as binned. Defaults to FALSE.
#'
#' @return A scMaSigProClass object containing the inputted counts, cell data, and additional parameters.
#'
#'
#' @export

# Create scmp
create.scmp <- function(counts,
                        cell_data,
                        feature_data,
                        bin_counts = NULL,
                        bin_cell_data = NULL,
                        pseudotime_colname,
                        path_colname,
                        use_as_bin = FALSE) {
  # Validation Checks
  assert_that(ncol(counts) == nrow(cell_data),
    msg = paste("Number of cells in raw-counts and cell-level-metadata are different.")
  )

  assert_that(all(colnames(counts) == rownames(cell_data)),
    msg = paste("Rownames of raw-counts and cell-level-metadata are different.")
  )

  if (!is.null(bin_counts) | !is.null(bin_cell_data)) {
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
  scmpObj <- new("scMaSigProClass",
    sparse = sparse_tmp,
    dense = SingleCellExperiment(assays = list(bulk.counts = matrix(0, nrow = 0, ncol = 0)))
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
    scmpObj@dense <- sparse_tmp
    sparse_tmp <- NULL

    # Update the slots
    scmpObj@param@bin_pseudotime_colname <- pseudotime_colname
  }

  # Update the slots
  scmpObj@param@pseudotime_colname <- pseudotime_colname
  scmpObj@param@path_colname <- path_colname


  return(scmpObj)
}
