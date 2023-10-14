#' Annotate SingleCellExperiment object with pseudotime and path information.
#'
#' This function allows you to annotate a SingleCellExperiment object with pseudotime
#' and path information in its cell metadata.
#'
#' @param sce A SingleCellExperiment object to be annotated.
#' @param pseudotime_colname The name of the pseudotime column to be added or updated.
#' @param path_prefix The prefix for the path column name (default is "Path").
#' @param root_label The label for the root of the path (default is "root").
#' @param path_colname The name of the path column to be added or updated.
#' @param existing_pseudotime_colname The name of an existing pseudotime column to be replaced (if not NULL).
#' @param existing_path_colname The name of an existing path column to be replaced (if not NULL).
#' @param overwrite_labels Logical, should existing column names be overwritten if they exist? (default is TRUE).
#' @param verbose Logical, should progress messages be printed? (default is TRUE).
#'
#' @return A SingleCellExperiment object with updated cell metadata.
#'
#' @examples
#' # Annotate a SingleCellExperiment object with pseudotime and path information
#' \dontrun{
#' annotated_sce <- annotate_sce(sce, pseudotime_colname = "Pseudotime", path_colname = "Path")
#' }
#'
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#' @importFrom SingleCellExperiment colData
#' @importFrom assertthat assert_that
#'
#' @export
annotate_sce <- function(sce,
                         pseudotime_colname = "Pseudotime",
                         path_prefix = "Path",
                         root_label = "root",
                         path_colname = path_prefix,
                         existing_pseudotime_colname = NULL,
                         existing_path_colname = NULL,
                         overwrite_labels = TRUE,
                         verbose = TRUE) {
  # Overwite the columns
  if (overwrite_labels) {
    assert_that(
      all(!is.null(existing_pseudotime_colname) & !is.null(existing_path_colname)),
      msg = "If 'overwrite_labels' is TRUE, 'existing_pseudotime_colname' and 'existing_path_colname', cannot be NULL"
    )
    
    # Extract the cell metadata
    cell.meta <- as.data.frame(colData(sce))
    
    # Check columns
    assert_that(
        all(existing_pseudotime_colname %in% colnames(cell.meta)),
        msg = "'existing_pseudotime_colname', doesn't exist in cell.metadata"
    )
    # Check columns
    assert_that(
        all(existing_path_colname %in% colnames(cell.meta)),
        msg = "'existing_path_colname', doesn't exist in cell.metadata"
    )
    
    # Override
    names(cell.meta)[names(cell.meta) == existing_pseudotime_colname] <- pseudotime_colname
    
    # Overwite the columns
    names(cell.meta)[names(cell.meta) == existing_path_colname] <- path_colname

    # Update cell dataset with the updated cell metadata
    colData(sce) <- DataFrame(cell.meta)

    # Return
    return(sce)
  } else {
    # Return
    return(sce)
  }
}
