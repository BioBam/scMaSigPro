#' Create scMaSigProClass Object
#'
#' This function initializes a scMaSigProClass object with the given counts, cell data, and other optional parameters.
#' It is designed for single-cell RNA-seq data analysis.
#'
#' @param counts A matrix containing the expression counts.
#' @param cell_data A data frame containing the cell metadata.
#' @param feature_data A data frame containing the feature metadata. (Currently unused)
#' @param bin_counts A matrix containing the binned counts. (Currently unused)
#' @param bin_cell_data A data frame containing the binned cell data. (Currently unused)
#' @param pseudotime_colname A character string specifying the column name for pseudotime in the cell metadata.
#' @param path_colname A character string specifying the column name for the path in the cell metadata.
#' @param use_as_bin A logical flag indicating whether to use the counts and cell data as bin. Defaults to FALSE.
#'
#' @return A scMaSigProClass object containing the inputted counts, cell data, and additional parameters.
#'
#' @examples
#' \dontrun{
#' # Assuming you have 'counts_matrix', 'cell_data_df', 'pseudotime_col', and 'path_col'
#' scmp_object <- create_scmp(counts_matrix, cell_data_df,
#'   feature_data = NULL,
#'   bin_counts = NULL, bin_cell_data = NULL,
#'   pseudotime_col, path_col, use_as_bin = FALSE
#' )
#' }
#'
#' @export
#'

# Create scmp
create_scmpObj <- function(counts, cell_data, feature_data,
                           bin_counts, bin_cell_data,
                           pseudotime_colname,
                           path_colname,
                           use_as_bin = FALSE) {
  # Create Single-Cell Experiment Object
  sce_tmp <- SingleCellExperiment(
    list(counts = counts),
    colData = DataFrame(cell_data)
  )

  # Initate scMaSigPro
  scmpObj <- new("scMaSigProClass",
    sce = sce_tmp,
    compress.sce = SingleCellExperiment(assays = list(bulk.counts = matrix(0, nrow = 0, ncol = 0)))
  )
  sce_tmp <- NULL

  # Use as bin
  if (use_as_bin) {
    sce_tmp <- SingleCellExperiment(
      list(bulk.counts = counts),
      colData = DataFrame(cell_data)
    )

    # Transfer Data
    scmpObj@compress.sce <- sce_tmp
    sce_tmp <- NULL

    # Update the slots
    scmpObj@addParams@bin_pseudotime_colname <- pseudotime_colname
  }

  # Update the slots
  scmpObj@addParams@pseudotime_colname <- pseudotime_colname
  scmpObj@addParams@path_colname <- path_colname


  return(scmpObj)
}
