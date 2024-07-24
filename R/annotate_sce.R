#' @title Annotate `SingleCellExperiment` class object with pseudotime and path
#' information.
#'
#' @description
#' `annotate_sce()` annotates a `SingleCellExperiment` class object with pseudotime
#' and path information in its `cell.metadata` generated using `colData` from the
#' \pkg{SingleCellExperiment} package.
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom assertthat assert_that
#'
#' @param sce A `SingleCellExperiment` object to be annotated.
#' @param ptime_col A character string representing the column name
#' for inferred Pseudotime values. (Default is "Pseudotime")
#' @param path_prefix Prefix used to annotate the paths. (Default is "Path").
#' @param root_label Label used to annotate root cells. (Default is "root").
#' @param path_col A character string representing the column name for branching
#' path assignment. (Default is `path_prefix`)
#' @param labels_exist Logical indicating whether if the existing column names
#' should be overwritten. (Default is TRUE)
#' @param exist_ptime_col The name of an existing pseudotime column
#' to be replaced (if not NULL).
#' @param exist_path_col The name of an existing path column to be replaced
#' (if not NULL).
#' @param verbose Print detailed output in the console. (Default is TRUE)
#'
#' @return A \code{\link{SingleCellExperiment}} object with updated cell metadata.
#'
#' @seealso \pkg{SingleCellExperiment} package.
#' @author Priyansh Srivastava \email{spriyansh29@@gmail.com}
#' @keywords internal
#'
annotate_sce <- function(sce,
                         ptime_col = "Pseudotime",
                         path_prefix = "Path",
                         root_label = "root",
                         path_col = path_prefix,
                         exist_ptime_col = NULL,
                         exist_path_col = NULL,
                         labels_exist = FALSE,
                         verbose = TRUE) {
  # Overwite the columns
  if (labels_exist) {
    assertthat::assert_that(
      all(!is.null(exist_ptime_col) & !is.null(exist_path_col)),
      msg = paste(
        "Requested to set 'path_col' as", path_col,
        "with 'labels_exist' as TRUE. Please supply
                  'exist_ptime_col' and 'exist_path_col'"
      )
    )

    if (verbose) {
      message(paste0(
        "Overwritting columns in cell level metadata, '", exist_path_col, "'
        is replaced by '", path_col,
        "' and '", exist_ptime_col, "' is replaced by '", ptime_col, "'."
      ))
    }

    # Extract the cell metadata
    cell.meta <- as.data.frame(colData(sce))

    # Check columns
    assertthat::assert_that(
      all(exist_ptime_col %in% colnames(cell.meta)),
      msg = paste("'", exist_ptime_col, "', doesn't exist in colData.")
    )
    # Check columns
    assertthat::assert_that(
      all(exist_path_col %in% colnames(cell.meta)),
      msg = paste("'", exist_path_col, "', doesn't exist in colData")
    )

    # Override
    names(cell.meta)[names(cell.meta) == exist_ptime_col] <- ptime_col

    # Overwite the columns
    names(cell.meta)[names(cell.meta) == exist_path_col] <- path_col

    # Update cell dataset with the updated cell metadata
    sce@colData <- DataFrame(cell.meta)

    # Return
    return(sce)
  } else {
    if (verbose) {
      message(paste("Skipping overwritting of colnames in cell level metadata"))
    }
    # Return
    return(sce)
  }
}
