#' Create scMaSigProClass Object
#'
#' This function initializes a scMaSigProClass object with the given counts,
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
create_scmpObj <- function(counts, 
                           cell_data,
                           feature_data,
                           bin_counts = NULL,
                           bin_cell_data = NULL,
                           pseudotime_colname,
                           path_colname,
                           use_as_bin = FALSE) {
    
    # Validation Checks
    assert_that(ncol(counts) == nrow(cell_data),
                msg = paste("Number of cells in raw-counts and cell-level-metadata are different."))
    
    assert_that(all(colnames(counts) == rownames(cell_data)),
                msg = paste("Rownames of raw-counts and cell-level-metadata are different."))
    
    if(!is.null(bin_counts) | !is.null(bin_cell_data)){
        assert_that(nrow(bin_counts) == nrow(bin_cell_data),
                    msg = paste("Number of cells in bin_counts and bin_cell_data are different."))
        assert_that(all(rownames(bin_counts) == rownames(bin_cell_data)),
                    msg = paste("Rownames of bin_counts and bin_cell_data are different."))
    }
    
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
